"""
Gait analysis pipeline for the Synced-IMU-Data dataset.

Processes synced subjects across the four reported conditions (UL, L2, L3,
L7) using:
  - the adaptive-threshold gait event detector
  - the ZUPT stride length estimator
  - GaitRite ground truth for accuracy validation
  - GR Pulse time synchronization for IMU-GaitRite alignment

Outputs (in {base_dir}/outputs/{date}/):
  - IMU_Synced_results.json                (per-subject processing summary)
  - IMU_Synced_accuracy_summary.csv        (detection accuracy by condition)
  - IMU_Synced_per_subject_detection.csv   (per-subject detection rates)
  - IMU_Synced_spatiotemporal_summary.csv  (spatiotemporal by condition)
  - IMU_Synced_stride_length_summary.csv   (stride length by condition)
  - stride_level_pairs.csv                 (per-stride ZUPT-GaitRite pairs)
  - alignment_qc.csv                       (walk-level alignment decisions)

Subjects in every output are labelled by their de-identified study code
(ID_NNN), derived from the raw data folder name by subject_code().

Usage:
    python run_pipeline.py --base-dir /path/to/Synced-IMU-Data
    python run_pipeline.py --subjects <folder name> --conditions UL L2
"""

import argparse
import csv
import json
import os
import re
import sys
from datetime import date
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from comparison import match_events
from config import (
    BANDPASS_HIGH,
    BANDPASS_LOW,
    CONDITIONS,
    EXCLUDE_STRIDE_PAIRS,
    EXCLUDE_WALKS,
    FILTER_ORDER,
    GR_SL_MAX_M,
    GR_SL_MIN_M,
    GR_ST_MAX_S,
    GR_ST_MIN_S,
    INVERT_LEFT_GYRO,
    INVERT_RIGHT_GYRO,
    LEFT_ACCEL_COLS,
    LEFT_GYRO_COLS,
    LEFT_GYRO_Z_COL,
    RIGHT_ACCEL_COLS,
    RIGHT_GYRO_COLS,
    RIGHT_GYRO_Z_COL,
    SAMPLING_FREQ,
    subject_code,
)
from gait_events import AdaptiveThresholdDetector
from stride_length import (
    estimate_stride_length_zupt,
)
from time_sync import (
    find_optimal_alignment,
    load_gr_pulse,
    split_gaitrite_data,
)

BASE_DIR = os.environ.get("IMU_SYNCED_BASE_DIR")

MATCH_TOL_S = 0.150
OUTLIER_FRACTION = 0.10

EVENT_KEYS = ["left_hs", "right_hs", "left_to", "right_to"]


def imu_path(subject, cond_code):
    folder = CONDITIONS[cond_code][0]
    return os.path.join(
        BASE_DIR, "IMU Data", subject, folder, f"{subject}_{cond_code}.csv"
    )


def gr_pulse_path(subject, cond_code):
    base = os.path.join(BASE_DIR, "Time Sychng Information Data", subject)
    # Try short code first (e.g. _UL_), then folder name (e.g. _Unloaded_)
    for label in [cond_code, CONDITIONS[cond_code][0]]:
        path = os.path.join(base, f"{subject}_{label}_GR_Pulse.csv")
        if os.path.exists(path):
            return path
    return os.path.join(base, f"{subject}_{cond_code}_GR_Pulse.csv")


def gaitrite_dir(subject, cond_code):
    folder = CONDITIONS[cond_code][0]
    root = os.path.join(BASE_DIR, "GaitRite Data", subject)
    # Historical exports use both "Loaded 2" and "Loaded2" folder names.
    for candidate in [folder, folder.replace(" ", "")]:
        path = os.path.join(root, candidate)
        if os.path.isdir(path):
            return path
    return os.path.join(root, folder)


def load_imu_csv(filepath):
    """Load and bandpass-filter an IMU CSV.

    Returns dict with time, filtered gyro-z signals, and raw 3-axis
    accel/gyro for ZUPT.
    """
    df = pd.read_csv(filepath, skiprows=3, header=None)

    time_us = df.iloc[:, 0].values
    time_s = (time_us - time_us[0]) / 1e6

    # Derive observed fs from timestamps and design the filter against it
    # rather than the nominal SAMPLING_FREQ. APDM Opal nominally streams at
    # 256 Hz but clock drift and buffer underruns can shift the effective
    # rate; designing the bandpass against the wrong Nyquist would shift
    # the cutoff in the wrong direction.
    fs_observed = SAMPLING_FREQ
    if len(time_s) > 1:
        dt = np.median(np.diff(time_s))
        if dt > 0:
            fs_observed = 1.0 / dt
            if abs(fs_observed - SAMPLING_FREQ) / SAMPLING_FREQ > 0.02:
                print(
                    f"    WARNING: observed fs={fs_observed:.1f} Hz "
                    f"differs from SAMPLING_FREQ={SAMPLING_FREQ} Hz by "
                    f">2% in {os.path.basename(filepath)}"
                )

    right_gyro_z = df.iloc[:, RIGHT_GYRO_Z_COL].values
    left_gyro_z = df.iloc[:, LEFT_GYRO_Z_COL].values

    # Sign convention for HS detection only. Do NOT propagate to the
    # 3-axis arrays loaded below: a single-axis flip would create a
    # left-handed frame and corrupt Mahony AHRS orientation tracking.
    if INVERT_LEFT_GYRO:
        left_gyro_z = -left_gyro_z
    if INVERT_RIGHT_GYRO:
        right_gyro_z = -right_gyro_z

    nyq = fs_observed / 2
    b, a = butter(FILTER_ORDER, [BANDPASS_LOW / nyq, BANDPASS_HIGH / nyq], btype="band")
    right_filt = filtfilt(b, a, right_gyro_z)
    left_filt = filtfilt(b, a, left_gyro_z)

    result = {
        "time": time_s,
        "right_filt": right_filt,
        "left_filt": left_filt,
        "fs": fs_observed,
    }

    ncols = df.shape[1]
    if max(RIGHT_ACCEL_COLS + RIGHT_GYRO_COLS) < ncols:
        result["right_accel_3ax"] = df.iloc[:, RIGHT_ACCEL_COLS].values
        result["right_gyro_3ax"] = df.iloc[:, RIGHT_GYRO_COLS].values
    if max(LEFT_ACCEL_COLS + LEFT_GYRO_COLS) < ncols:
        result["left_accel_3ax"] = df.iloc[:, LEFT_ACCEL_COLS].values
        result["left_gyro_3ax"] = df.iloc[:, LEFT_GYRO_COLS].values

    return result


def load_trial_data(subject, cond_code, imu_data):
    """Load GR Pulse walks and GaitRite trial files for one subject/condition.

    Returns list of trial dicts, each containing trial_num, gr_walks
    (OUT/BACK), gr_spatial, tdms_out, tdms_back.
    """
    pulse_file = gr_pulse_path(subject, cond_code)
    if not os.path.exists(pulse_file):
        raise FileNotFoundError(f"GR Pulse not found: {pulse_file}")

    tdms_walks = load_gr_pulse(pulse_file, imu_data=imu_data)

    gr_dir = gaitrite_dir(subject, cond_code)
    if not os.path.isdir(gr_dir):
        raise FileNotFoundError(f"GaitRite dir not found: {gr_dir}")

    xlsx_files = sorted(glob(os.path.join(gr_dir, "*.xlsx")))
    trials = []
    for filepath in xlsx_files:
        filename = Path(filepath).name
        match = re.search(r"Trial_(\d+)", filename)
        if not match:
            continue
        trial_num = int(match.group(1))
        try:
            gr_walks = split_gaitrite_data(filepath)
        except Exception as e:
            print(f"    Warning: {filename}: {e}")
            continue

        # Each trial maps to 2 TDMS walks: out and back
        tdms_idx = (trial_num - 1) * 2
        if tdms_idx + 1 >= len(tdms_walks):
            continue

        try:
            gr_spatial = extract_gaitrite_spatial(filepath)
        except Exception as e:
            print(f"    Warning: spatial extraction failed for {filename}: {e}")
            gr_spatial = {}

        trials.append(
            {
                "trial_num": trial_num,
                "gr_walks": gr_walks,
                "gr_spatial": gr_spatial,
                "tdms_out": tdms_walks[tdms_idx],
                "tdms_back": tdms_walks[tdms_idx + 1],
            }
        )

    return trials


def extract_gaitrite_spatial(filepath):
    """Extract per-walk stride data from a GaitRite xlsx file.

    GaitRite reports per-footfall data. Each trial contains an OUT walk
    followed by a BACK walk. The first two footfalls per walk lack stride
    data (a stride needs 2 contacts).

    Returns dict with keys 'OUT' and 'BACK', each containing
    stride_length_cm, stride_time, and foot (0=left, 1=right) arrays.
    """
    df = pd.read_excel(filepath)
    df = df.dropna(subset=["FootFall Object #"])

    # Split into OUT and BACK where FootFall Object # resets to 1
    split_idx = None
    ffo = df["FootFall Object #"].values
    for i in range(1, len(ffo)):
        if ffo[i] == 1.0:
            split_idx = i
            break

    if split_idx is None:
        walks_df = {"OUT": df}
    else:
        walks_df = {
            "OUT": df.iloc[:split_idx].copy(),
            "BACK": df.iloc[split_idx:].copy(),
        }

    result = {}
    for walk_key, wdf in walks_df.items():
        valid = wdf[wdf["Stride Length"] > 0].copy()
        result[walk_key] = {
            "stride_length_cm": valid["Stride Length"].values,
            "stride_time": valid["Stride Time"].values,
            "foot": valid["Left/Right Foot"].values,
        }

    return result


def bidirectional_filter(imu_events, gr_events, fraction=OUTLIER_FRACTION):
    """Remove outlier events that have no nearby match in the other set.

    Tolerance uses a fixed floor of MATCH_TOL_S (150 ms). The previous
    `median_stride_interval * 0.10` (~50-100 ms) was tighter than the
    detector RMSE for irregular gait, silently dropping the dysrhythmic
    strides most relevant to validation in MCI / dual-task subjects.
    """
    if len(gr_events) < 2 or len(imu_events) == 0:
        return imu_events, gr_events
    tol = max(MATCH_TOL_S, np.median(np.diff(np.sort(gr_events))) * fraction)

    gr_keep = [
        g
        for g in gr_events
        if len(imu_events) > 0 and np.min(np.abs(imu_events - g)) <= tol
    ]
    gr_clean = np.array(gr_keep) if gr_keep else np.array([])
    if len(gr_clean) == 0:
        return np.array([]), np.array([])

    imu_keep = [i for i in imu_events if np.min(np.abs(gr_clean - i)) <= tol]
    return (np.array(imu_keep) if imu_keep else np.array([])), gr_clean


def compute_spatiotemporal(left_hs, right_hs, left_to, right_to):
    """Compute temporal gait parameters from detected events."""
    results = {}

    for side, hs in [("left", left_hs), ("right", right_hs)]:
        s = np.sort(hs)
        strides = np.diff(s) if len(s) >= 2 else np.array([])
        valid = (
            strides[(strides >= 0.4) & (strides <= 3.0)]
            if len(strides) > 0
            else np.array([])
        )
        results[f"{side}_stride_time"] = valid

    all_hs = sorted([(t, "L") for t in left_hs] + [(t, "R") for t in right_hs])
    steps = []
    for i in range(1, len(all_hs)):
        if all_hs[i][1] != all_hs[i - 1][1]:
            st = all_hs[i][0] - all_hs[i - 1][0]
            if 0.2 <= st <= 2.0:
                steps.append(st)
    results["step_time"] = np.array(steps)

    for side, hs, to in [("left", left_hs, left_to), ("right", right_hs, right_to)]:
        hs_s, to_s = np.sort(hs), np.sort(to)
        stance = []
        for h in hs_s:
            cands = to_s[to_s > h]
            if len(cands) > 0:
                dur = cands[0] - h
                if 0.1 <= dur <= 2.0:
                    stance.append(dur)
        results[f"{side}_stance_time"] = np.array(stance)

    for side, hs, to in [("left", left_hs, left_to), ("right", right_hs, right_to)]:
        hs_s, to_s = np.sort(hs), np.sort(to)
        swing = []
        for t_off in to_s:
            cands = hs_s[hs_s > t_off]
            if len(cands) > 0:
                dur = cands[0] - t_off
                if 0.1 <= dur <= 2.0:
                    swing.append(dur)
        results[f"{side}_swing_time"] = np.array(swing)

    return results


def summarize_array(arr):
    """Mean +/- SD summary."""
    if len(arr) == 0:
        return {"mean": None, "sd": None, "n": 0}
    return {
        "mean": round(float(np.mean(arr)), 3),
        "sd": round(float(np.std(arr, ddof=1)), 3) if len(arr) > 1 else 0.0,
        "n": int(len(arr)),
    }


def discover_subjects():
    """Find all subject folders under IMU Data."""
    imu_dir = os.path.join(BASE_DIR, "IMU Data")
    subjects = []
    for entry in sorted(os.listdir(imu_dir)):
        full = os.path.join(imu_dir, entry)
        if os.path.isdir(full):
            subjects.append(entry)
    return subjects


def main():
    global BASE_DIR
    parser = argparse.ArgumentParser(
        description="Gait analysis pipeline for the Synced-IMU-Data dataset"
    )
    parser.add_argument(
        "--base-dir",
        default=BASE_DIR,
        help="Synced-IMU-Data root (or set IMU_SYNCED_BASE_DIR)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: BASE_DIR/outputs/YYYY-MM-DD)",
    )
    parser.add_argument(
        "--subjects",
        nargs="+",
        default=None,
        help="Subject IDs to process (default: all)",
    )
    parser.add_argument(
        "--conditions",
        nargs="+",
        default=None,
        help="Condition codes to process (default: all)",
    )
    args = parser.parse_args()

    if not args.base_dir:
        parser.error("provide --base-dir or set IMU_SYNCED_BASE_DIR")
    BASE_DIR = os.path.abspath(os.path.expanduser(args.base_dir))

    # The adaptive detector derives its thresholds from each walk's own
    # gyro amplitude, so it stays reliable for the slow, low-amplitude
    # shuffling gait seen in some older-adult and dual-task trials where
    # a fixed-threshold detector would fail silently. It supplies both
    # the reported gait events and the IMU-side anchor for sync alignment.
    detector = AdaptiveThresholdDetector()
    subjects = discover_subjects()

    if args.subjects:
        requested = set(args.subjects)
        subjects = [s for s in subjects if s in requested]
        missing = requested - set(subjects)
        if missing:
            print(f"WARNING: subjects not found: {missing}")

    cond_list = ["UL", "L2", "L3", "L7"]
    if args.conditions:
        cond_list = [c for c in args.conditions if c in CONDITIONS]

    print(f"Base directory: {BASE_DIR}")
    print(f"Subjects: {len(subjects)}")
    print(f"Conditions: {cond_list}")

    # accuracy: cond -> event_key -> list of error_ms
    acc_errors = {}
    acc_matched = {}
    acc_ref = {}
    acc_fp = {}

    # spatiotemporal: cond -> source (IMU / GaitRite) -> param -> [values]
    st_accum = {}

    # stride length: cond -> side -> [values]
    sl_accum = {}

    all_results = {}  # subject_cond -> summary dict
    per_subj_event_rows = []
    stride_level_rows = []
    alignment_qc_rows = []

    for cond_code in cond_list:
        cond_folder, cond_label = CONDITIONS[cond_code]
        print(f"\nCONDITION: {cond_folder} ({cond_code})")

        acc_errors[cond_code] = {ek: [] for ek in EVENT_KEYS}
        acc_matched[cond_code] = {ek: 0 for ek in EVENT_KEYS}
        acc_ref[cond_code] = {ek: 0 for ek in EVENT_KEYS}
        acc_fp[cond_code] = {ek: 0 for ek in EVENT_KEYS}

        st_params = [
            "left_stride_time",
            "right_stride_time",
            "step_time",
            "left_stance_time",
            "right_stance_time",
            "left_swing_time",
            "right_swing_time",
        ]
        st_accum[cond_code] = {
            src: {p: [] for p in st_params} for src in ["IMU", "GaitRite"]
        }
        sl_accum[cond_code] = {"left_zupt": [], "right_zupt": []}

        subjects_ok = 0
        total_walks = 0

        for subject_id in subjects:
            imu_file = imu_path(subject_id, cond_code)
            if not os.path.exists(imu_file):
                continue

            print(f"\n  {subject_id} ({cond_code})")

            try:
                imu = load_imu_csv(imu_file)
            except Exception as e:
                print(f"    ERROR loading IMU: {e}")
                continue

            try:
                trials = load_trial_data(subject_id, cond_code, imu)
            except Exception as e:
                print(f"    ERROR loading trials: {e}")
                continue

            if len(trials) == 0:
                print("    No valid trials found")
                continue

            subjects_ok += 1
            subject_label = subject_code(subject_id)
            subj_matched = {ek: 0 for ek in EVENT_KEYS}
            subj_ref = {ek: 0 for ek in EVENT_KEYS}
            time_vec = imu["time"]
            left_filt = imu["left_filt"]
            right_filt = imu["right_filt"]
            fs = imu.get("fs", SAMPLING_FREQ)

            alignment_cache = {}

            def get_walk_alignment(trial, walk_key, tdms_walk):
                """Return one cached, quality-checked GR-to-IMU alignment."""
                cache_key = (trial["trial_num"], walk_key)
                if cache_key in alignment_cache:
                    return alignment_cache[cache_key]

                gr_events = trial.get("gr_walks", {}).get(walk_key, {})
                gr_left = np.asarray(gr_events.get("left_hs", []))
                gr_right = np.asarray(gr_events.get("right_hs", []))
                imu_start = tdms_walk.get("imu_start")
                imu_end = tdms_walk.get("imu_end")
                result = {
                    "valid": False,
                    "offset_s": np.nan,
                    "swap_suggested": False,
                    "mean_abs_error_ms": np.nan,
                    "n_matched": 0,
                    "n_expected": int(len(gr_left) + len(gr_right)),
                    "match_ratio": 0.0,
                    "imu_left_hs": 0,
                    "imu_right_hs": 0,
                    "reason": "",
                }

                if imu_start is None or imu_end is None:
                    result["reason"] = "Missing GR Pulse walk boundary"
                elif result["n_expected"] == 0:
                    result["reason"] = "No GaitRite heel strikes"
                else:
                    mask_w = (time_vec >= imu_start - 2.0) & (time_vec <= imu_end + 1.5)
                    t_local = time_vec[mask_w]
                    if len(t_local) < fs:
                        result["reason"] = "IMU walk window too short"
                    else:
                        try:
                            align_lhs = detector.detect_heel_strikes(
                                t_local, left_filt[mask_w], fs=fs
                            )
                            align_rhs = detector.detect_heel_strikes(
                                t_local, right_filt[mask_w], fs=fs
                            )
                            result["imu_left_hs"] = int(len(align_lhs))
                            result["imu_right_hs"] = int(len(align_rhs))
                            offset, swap, stats = find_optimal_alignment(
                                align_lhs, align_rhs, gr_left, gr_right
                            )
                            n_matched = int(stats.get("n_matched", 0))
                            n_expected = int(
                                stats.get("n_expected", result["n_expected"])
                            )
                            mean_error = stats.get("mean_abs_error_ms", np.nan)
                            match_ratio = n_matched / n_expected if n_expected else 0.0
                            valid = (
                                "error" not in stats
                                and np.isfinite(offset)
                                and np.isfinite(mean_error)
                                and n_matched > 0
                                and match_ratio >= 0.60
                            )
                            result.update(
                                {
                                    "valid": bool(valid),
                                    "offset_s": float(offset),
                                    "swap_suggested": bool(swap),
                                    "mean_abs_error_ms": (
                                        float(mean_error)
                                        if np.isfinite(mean_error)
                                        else np.nan
                                    ),
                                    "n_matched": n_matched,
                                    "n_expected": n_expected,
                                    "match_ratio": float(match_ratio),
                                    "reason": (
                                        ""
                                        if valid
                                        else stats.get(
                                            "error", "Alignment quality check failed"
                                        )
                                    ),
                                }
                            )
                        except (ValueError, IndexError, RuntimeError) as exc:
                            result["reason"] = f"Alignment error: {exc}"

                alignment_cache[cache_key] = result
                alignment_qc_rows.append(
                    {
                        "subject": subject_label,
                        "condition": cond_code,
                        "trial_num": trial["trial_num"],
                        "walk": walk_key,
                        **result,
                    }
                )
                if not result["valid"]:
                    print(
                        f"    SYNC QC: excluding {subject_id} {cond_code} "
                        f"trial {trial['trial_num']} {walk_key} from paired "
                        f"validation ({result['reason']})"
                    )
                return result

            # Run ZUPT once per foot for this subject/condition. Keep both
            # stride lengths and their timestamps so we can filter to
            # GaitRite on-mat windows later.
            zupt_raw = {}
            side_config = {
                "left": {"ipsi_filt": left_filt, "contra_filt": right_filt},
                "right": {"ipsi_filt": right_filt, "contra_filt": left_filt},
            }
            for side in ["left", "right"]:
                accel_key = f"{side}_accel_3ax"
                gyro_key = f"{side}_gyro_3ax"
                if accel_key in imu and gyro_key in imu:
                    try:
                        zupt_out = estimate_stride_length_zupt(
                            imu[accel_key],
                            imu[gyro_key],
                            time_vec,
                            fs,
                            ipsi_gyro_z_filt=side_config[side]["ipsi_filt"],
                            contra_gyro_z_filt=side_config[side]["contra_filt"],
                            stance_method="contralateral",
                            bidirectional=True,
                        )
                        zupt_raw[side] = zupt_out
                    except Exception as e:
                        print(f"    [ZUPT {side}] Error: {e}")
                        zupt_raw[side] = None
                else:
                    zupt_raw[side] = None

            # Collect GaitRite on-mat time windows from GR Pulse walks. Each
            # walk boundary (imu_start to imu_end) defines when the person
            # is walking across the mat; strides outside these windows are
            # turnarounds, rest, or non-GR segments.
            on_mat_windows = []
            for trial in trials:
                for tdms_walk in [trial["tdms_out"], trial["tdms_back"]]:
                    w_start = tdms_walk.get("imu_start")
                    w_end = tdms_walk.get("imu_end")
                    if w_start is not None and w_end is not None:
                        on_mat_windows.append((w_start, w_end))

            # Filter ZUPT strides to on-mat windows only
            zupt_results = {}
            zupt_results_with_times = {}
            for side in ["left", "right"]:
                zr = zupt_raw.get(side)
                if zr is None or zr["n_strides"] == 0:
                    zupt_results[side] = np.array([])
                    zupt_results_with_times[side] = {
                        "sl": np.array([]),
                        "times": np.array([]),
                    }
                    continue

                sl = zr["stride_lengths"]
                # ZUPT reports stride_midpoints anchored at stance-segment
                # start (HS-equivalent), matching GaitRite's
                # (HS_i + HS_{i+1})/2 stride midpoint.
                stride_times_abs = zr["stride_midpoints"]
                st_dur = zr["stride_times"]

                # Keep strides whose full HS-to-HS window falls inside an
                # on-mat window. Filtering by midpoint alone admits strides
                # that straddle the mat boundary (e.g. a heel-strike before
                # the foot enters the mat), which GR cannot record and so
                # produces a spurious pair against an unrelated GR stride.
                hs_start = stride_times_abs - st_dur / 2.0
                hs_end = stride_times_abs + st_dur / 2.0
                keep = np.zeros(len(sl), dtype=bool)
                for w_start, w_end in on_mat_windows:
                    keep |= (hs_start >= w_start) & (hs_end <= w_end)

                mask = keep & np.isfinite(sl)
                valid = sl[mask]
                valid_times = stride_times_abs[mask]
                valid_durations = st_dur[mask]
                zupt_results[side] = valid
                zupt_results_with_times[side] = {
                    "sl": valid,
                    "times": valid_times,
                    "st": valid_durations,
                }

            lz = zupt_results.get("left", np.array([]))
            rz = zupt_results.get("right", np.array([]))
            n_zupt = len(lz) + len(rz)
            if n_zupt > 0:
                all_z = np.concatenate([lz, rz])
                print(
                    f"    ZUPT: {n_zupt} valid on-mat strides "
                    f"(L={len(lz)}, R={len(rz)}), "
                    f"mean={np.mean(all_z):.3f} m"
                )
            else:
                print("    ZUPT: 0 valid on-mat strides")

            # Stride-level pairing: match ZUPT strides to GR strides by
            # temporal proximity.
            if (subject_label, cond_code) in EXCLUDE_STRIDE_PAIRS:
                print("    EXCLUDED from stride pairing (known sync failure)")

            for trial in trials:
                if (subject_label, cond_code) in EXCLUDE_STRIDE_PAIRS:
                    break
                for gr_key, tdms_walk_key in [
                    ("OUT", "tdms_out"),
                    ("BACK", "tdms_back"),
                ]:
                    # Walk-level exclusion (data-quality issues documented
                    # in EXCLUDE_WALKS in config.py)
                    if (
                        subject_label,
                        cond_code,
                        trial["trial_num"],
                        gr_key,
                    ) in EXCLUDE_WALKS:
                        print(
                            f"    EXCLUDED walk {trial['trial_num']} {gr_key} "
                            f"(documented data-quality issue)"
                        )
                        continue
                    gr_sp = trial.get("gr_spatial", {}).get(gr_key, {})
                    gr_sl_cm = gr_sp.get("stride_length_cm", np.array([]))
                    gr_st_arr = gr_sp.get("stride_time", np.array([]))
                    gr_foot_arr = gr_sp.get("foot", np.array([]))
                    tdms_w = trial[tdms_walk_key]
                    if len(gr_sl_cm) == 0 or tdms_w.get("imu_start") is None:
                        continue

                    imu_s = tdms_w["imu_start"]
                    imu_e = tdms_w["imu_end"]

                    # GR events are in GaitRite time, not IMU time. A walk
                    # contributes only when heel-strike alignment passes QC.
                    gr_events_raw = trial.get("gr_walks", {}).get(gr_key, {})
                    gr_raw_hs = {
                        "left_hs": np.array(gr_events_raw.get("left_hs", [])),
                        "right_hs": np.array(gr_events_raw.get("right_hs", [])),
                    }
                    if (
                        len(gr_raw_hs["left_hs"]) == 0
                        and len(gr_raw_hs["right_hs"]) == 0
                    ):
                        continue

                    alignment = get_walk_alignment(trial, gr_key, tdms_w)
                    if not alignment["valid"]:
                        continue
                    walk_offset = alignment["offset_s"]

                    for side in ["left", "right"]:
                        zdata = zupt_results_with_times[side]
                        if len(zdata["sl"]) == 0:
                            continue

                        in_walk = (zdata["times"] >= imu_s) & (zdata["times"] <= imu_e)
                        z_sl_walk = zdata["sl"][in_walk]
                        z_t_walk = zdata["times"][in_walk]
                        z_st_walk = zdata["st"][in_walk]
                        if len(z_sl_walk) == 0:
                            continue

                        foot_val = 0 if side == "left" else 1
                        foot_mask = gr_foot_arr == foot_val
                        gr_sl_foot = gr_sl_cm[foot_mask] / 100.0  # cm -> m
                        gr_st_foot = (
                            gr_st_arr[foot_mask] if len(gr_st_arr) > 0 else np.array([])
                        )

                        if len(gr_sl_foot) == 0:
                            continue

                        # GR stride midpoint times from heel strikes,
                        # converted to IMU time using the alignment offset.
                        hs_key = f"{side}_hs"
                        gr_hs_raw = np.sort(np.array(gr_events_raw.get(hs_key, [])))
                        gr_hs = gr_hs_raw + walk_offset
                        n_gr_strides = len(gr_sl_foot)

                        if len(gr_hs) >= 2 and n_gr_strides > 0:
                            # GR stride i has midpoint at (hs[i] + hs[i+1]) / 2
                            n_hs_strides = len(gr_hs) - 1
                            n_use = min(n_gr_strides, n_hs_strides)
                            gr_midpoints = np.array(
                                [(gr_hs[i] + gr_hs[i + 1]) / 2.0 for i in range(n_use)]
                            )
                            gr_sl_timed = gr_sl_foot[:n_use]
                            gr_st_timed = (
                                gr_st_foot[:n_use]
                                if len(gr_st_foot) >= n_use
                                else np.array([])
                            )
                        else:
                            gr_midpoints = np.array([])
                            gr_sl_timed = np.array([])
                            gr_st_timed = np.array([])

                        if len(gr_midpoints) == 0:
                            continue

                        # GR reference plausibility: stride length AND stride
                        # time must fall within physiological bounds. The
                        # stride-time check rejects "double-stride" artifacts
                        # caused by GR missing an intermediate heel strike
                        # (typically light landings in slow / MCI walkers),
                        # which would otherwise be paired against a correctly
                        # detected IMU stride and inflate the apparent error.
                        plausible_sl = (gr_sl_timed >= GR_SL_MIN_M) & (
                            gr_sl_timed <= GR_SL_MAX_M
                        )
                        if (
                            len(gr_st_timed) == len(gr_sl_timed)
                            and len(gr_st_timed) > 0
                        ):
                            plausible_st = (gr_st_timed >= GR_ST_MIN_S) & (
                                gr_st_timed <= GR_ST_MAX_S
                            )
                            plausible = plausible_sl & plausible_st
                        else:
                            plausible = plausible_sl
                        gr_midpoints = gr_midpoints[plausible]
                        gr_sl_timed = gr_sl_timed[plausible]
                        gr_st_timed = (
                            gr_st_timed[plausible]
                            if len(gr_st_timed) > 0
                            else gr_st_timed
                        )
                        if len(gr_midpoints) == 0:
                            continue

                        # No first/last trim. Every stride that reaches this
                        # point is structurally complete:
                        #   - IMU: ZUPT requires two bounding stance phases
                        #   - GR: GaitRite reports stride_length only when both
                        #     heel strikes are captured on the mat
                        #   - On-mat filter: stride midpoint within [imu_s, imu_e]
                        #   - Mutual-nearest matching (below) rejects mis-pairs
                        z_sl_trim = z_sl_walk
                        z_t_trim = z_t_walk
                        z_st_trim = z_st_walk

                        # Mutual-nearest matching with cadence-scaled
                        # tolerance. A pair (gi, zi) is accepted only if gi
                        # is zi's nearest GR AND zi is gi's nearest IMU
                        # within tolerance, which avoids the cascading
                        # errors a fixed 0.5s tolerance produces after a
                        # single missed GR stride in slow L7 walkers.
                        if len(gr_midpoints) >= 2:
                            median_gr_dt = float(np.median(np.diff(gr_midpoints)))
                            match_tol = min(0.5, 0.4 * median_gr_dt)
                        else:
                            match_tol = 0.5

                        if len(z_t_trim) > 0 and len(gr_midpoints) > 0:
                            diff_matrix = np.abs(
                                gr_midpoints[:, None] - z_t_trim[None, :]
                            )
                            z_for_g = np.argmin(diff_matrix, axis=1)
                            g_for_z = np.argmin(diff_matrix, axis=0)
                        else:
                            z_for_g = np.array([], dtype=int)
                            g_for_z = np.array([], dtype=int)

                        for gi in range(len(gr_midpoints)):
                            if len(z_t_trim) == 0:
                                break
                            zi = int(z_for_g[gi])
                            if g_for_z[zi] != gi:
                                continue
                            diff_gz = float(diff_matrix[gi, zi])
                            if diff_gz > match_tol:
                                continue

                            # Plausibility filter: reject pairs whose ZUPT
                            # and GR stride length OR stride time disagree
                            # by >50%. This catches stance-detection
                            # failures where ZUPT misses a stance
                            # (integrates 2 strides, ratio ~2) or detects an
                            # extra stance (splits 1, ratio ~0.5). Such
                            # pairs are algorithmic failures, not
                            # measurement noise, and must not be reported
                            # as validation evidence.
                            zupt_st_val = float(z_st_trim[zi])
                            gr_sl_val = float(gr_sl_timed[gi])
                            gr_st_val = (
                                float(gr_st_timed[gi])
                                if gi < len(gr_st_timed)
                                else np.nan
                            )
                            sl_ratio = (
                                float(z_sl_trim[zi]) / gr_sl_val
                                if gr_sl_val > 0
                                else np.nan
                            )
                            st_ratio = (
                                zupt_st_val / gr_st_val
                                if (np.isfinite(gr_st_val) and gr_st_val > 0)
                                else np.nan
                            )
                            if not (0.5 <= sl_ratio <= 1.5):
                                continue
                            if np.isfinite(st_ratio) and not (0.5 <= st_ratio <= 1.5):
                                continue

                            zupt_speed = (
                                round(z_sl_trim[zi] / zupt_st_val, 4)
                                if zupt_st_val > 0
                                else np.nan
                            )
                            row = {
                                "subject": subject_label,
                                "condition": cond_code,
                                "condition_label": cond_folder,
                                "side": side,
                                "zupt_sl_m": round(float(z_sl_trim[zi]), 6),
                                "gr_sl_m": round(float(gr_sl_timed[gi]), 6),
                                "zupt_stride_time_s": round(zupt_st_val, 4),
                                "gr_stride_time_s": round(float(gr_st_timed[gi]), 4)
                                if gi < len(gr_st_timed)
                                else np.nan,
                                "zupt_speed_m_s": zupt_speed,
                                "zupt_time_s": round(float(z_t_trim[zi]), 4),
                                "gr_time_s": round(float(gr_midpoints[gi]), 4),
                                "time_diff_s": round(diff_gz, 4),
                                "trial_num": trial["trial_num"],
                                "walk": gr_key,
                            }
                            stride_level_rows.append(row)

            subject_walks = 0

            for trial in trials:
                for direction, tdms_walk, gr_key in [
                    ("OUT", trial["tdms_out"], "OUT"),
                    ("BACK", trial["tdms_back"], "BACK"),
                ]:
                    gr_data = trial["gr_walks"].get(gr_key)
                    if gr_data is None:
                        continue

                    gr_raw = {ek: np.array(gr_data.get(ek, [])) for ek in EVENT_KEYS}
                    if len(gr_raw["left_hs"]) == 0 and len(gr_raw["right_hs"]) == 0:
                        continue

                    imu_start = tdms_walk["imu_start"]
                    imu_end = tdms_walk["imu_end"]
                    mask = (time_vec >= imu_start - 2.0) & (time_vec <= imu_end + 1.5)
                    t_local = time_vec[mask]
                    left_local = left_filt[mask]
                    right_local = right_filt[mask]

                    if len(t_local) < fs:
                        continue

                    try:
                        lhs = detector.detect_heel_strikes(t_local, left_local, fs=fs)
                        rhs = detector.detect_heel_strikes(t_local, right_local, fs=fs)
                        lto = detector.detect_toe_offs(t_local, left_local, lhs, fs=fs)
                        rto = detector.detect_toe_offs(t_local, right_local, rhs, fs=fs)
                        events = {
                            "left_hs": lhs,
                            "right_hs": rhs,
                            "left_to": lto,
                            "right_to": rto,
                        }
                    except Exception as e:
                        print(f"    Detector error: {e}")
                        events = {ek: np.array([]) for ek in EVENT_KEYS}

                    alignment = get_walk_alignment(trial, gr_key, tdms_walk)
                    if not alignment["valid"]:
                        continue
                    best_offset = alignment["offset_s"]

                    gr_aligned = {ek: gr_raw[ek] + best_offset for ek in EVENT_KEYS}
                    total_walks += 1
                    subject_walks += 1

                    all_gr = np.concatenate(
                        [gr_aligned[ek] for ek in EVENT_KEYS if len(gr_aligned[ek]) > 0]
                    )
                    if len(all_gr) == 0:
                        continue
                    on_start = np.min(all_gr) - 0.3
                    on_end = np.max(all_gr) + 0.3

                    imu_filt = {}
                    gr_filt = {}
                    for ek in EVENT_KEYS:
                        imu_on = events[ek][
                            (events[ek] >= on_start) & (events[ek] <= on_end)
                        ]
                        imu_c, gr_c = bidirectional_filter(imu_on, gr_aligned[ek])
                        imu_filt[ek] = imu_c
                        gr_filt[ek] = gr_c

                    for ek in EVENT_KEYS:
                        matched_pairs, missed, false_pos = match_events(
                            imu_filt[ek], gr_filt[ek], MATCH_TOL_S
                        )
                        errs = [(det - ref) * 1000.0 for ref, det in matched_pairs]
                        acc_errors[cond_code][ek].extend(errs)
                        acc_matched[cond_code][ek] += len(matched_pairs)
                        acc_ref[cond_code][ek] += len(gr_filt[ek])
                        acc_fp[cond_code][ek] += len(false_pos)
                        subj_matched[ek] += len(matched_pairs)
                        subj_ref[ek] += len(gr_filt[ek])

                    imu_st = compute_spatiotemporal(
                        imu_filt["left_hs"],
                        imu_filt["right_hs"],
                        imu_filt["left_to"],
                        imu_filt["right_to"],
                    )
                    for key in st_params:
                        if key in imu_st and len(imu_st[key]) > 0:
                            st_accum[cond_code]["IMU"][key].extend(imu_st[key].tolist())

                    # GaitRite spatiotemporal from the same bidirectionally
                    # filtered reference events
                    gr_st = compute_spatiotemporal(
                        gr_filt["left_hs"],
                        gr_filt["right_hs"],
                        gr_filt["left_to"],
                        gr_filt["right_to"],
                    )
                    for key in st_params:
                        if key in gr_st and len(gr_st[key]) > 0:
                            st_accum[cond_code]["GaitRite"][key].extend(
                                gr_st[key].tolist()
                            )

            sl_accum[cond_code]["left_zupt"].extend(lz.tolist())
            sl_accum[cond_code]["right_zupt"].extend(rz.tolist())

            key = f"{subject_label}_{cond_code}"
            all_results[key] = {
                "subject": subject_label,
                "condition": cond_code,
                "n_trials": len(trials),
                "n_walks": subject_walks,
                "zupt_n_left": len(lz),
                "zupt_n_right": len(rz),
                "zupt_mean": round(float(np.mean(np.concatenate([lz, rz]))), 3)
                if n_zupt > 0
                else None,
            }

            print(f"    Trials: {len(trials)}, Walks: {subject_walks}")

            for event_label, keys in [
                ("HS", ["left_hs", "right_hs"]),
                ("TO", ["left_to", "right_to"]),
                ("Left HS", ["left_hs"]),
                ("Right HS", ["right_hs"]),
                ("Left TO", ["left_to"]),
                ("Right TO", ["right_to"]),
            ]:
                nm = sum(subj_matched[k] for k in keys)
                nr = sum(subj_ref[k] for k in keys)
                det_pct = 100.0 * nm / nr if nr > 0 else 0.0
                per_subj_event_rows.append(
                    {
                        "subject": subject_label,
                        "condition": cond_folder,
                        "event": event_label,
                        "n_matched": nm,
                        "n_ref": nr,
                        "detection_pct": round(det_pct, 2),
                    }
                )

        print(f"\n  Condition {cond_code}: {subjects_ok} subjects, {total_walks} walks")

    today = date.today().isoformat()
    output_dir = (
        os.path.abspath(os.path.expanduser(args.output_dir))
        if args.output_dir
        else os.path.join(BASE_DIR, "outputs", today)
    )
    os.makedirs(output_dir, exist_ok=True)

    # Detection accuracy summary, per-leg and pooled
    print("\nDETECTION ACCURACY BY CONDITION")

    accuracy_rows = []
    leg_combos = [
        ("Left HS", ["left_hs"]),
        ("Right HS", ["right_hs"]),
        ("HS", ["left_hs", "right_hs"]),
        ("Left TO", ["left_to"]),
        ("Right TO", ["right_to"]),
        ("TO", ["left_to", "right_to"]),
    ]

    for cond_code in cond_list:
        cond_folder = CONDITIONS[cond_code][0]

        print(f"\n  {cond_folder}")
        print(
            f"  {'Event':<10} {'Det%':>6} {'N_match':>8} {'N_ref':>7} "
            f"{'FP':>5} {'Bias':>7} {'SD':>7} {'RMSE':>7} {'MAE':>7}"
        )
        for label, keys in leg_combos:
            errs = []
            nm = nr = nf = 0
            for ek in keys:
                errs.extend(acc_errors[cond_code][ek])
                nm += acc_matched[cond_code][ek]
                nr += acc_ref[cond_code][ek]
                nf += acc_fp[cond_code][ek]

            errs = np.array(errs)
            if len(errs) > 0:
                bias = float(np.mean(errs))
                sd = float(np.std(errs, ddof=1))
                rmse = float(np.sqrt(np.mean(errs**2)))
                mae = float(np.mean(np.abs(errs)))
                det = 100.0 * nm / nr if nr > 0 else 0.0
            else:
                bias = sd = rmse = mae = det = 0.0

            print(
                f"  {label:<10} {det:>6.1f} {nm:>8} {nr:>7} "
                f"{nf:>5} {bias:>7.1f} {sd:>7.1f} {rmse:>7.1f} {mae:>7.1f}"
            )

            accuracy_rows.append(
                [
                    cond_folder,
                    label,
                    round(det, 1),
                    nm,
                    nr,
                    nf,
                    round(bias, 1),
                    round(sd, 1),
                    round(rmse, 1),
                    round(mae, 1),
                ]
            )

    acc_csv = os.path.join(output_dir, "IMU_Synced_accuracy_summary.csv")
    with open(acc_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "Condition",
                "Event",
                "Detection%",
                "N_matched",
                "N_ref",
                "N_false_pos",
                "Bias_ms",
                "SD_ms",
                "RMSE_ms",
                "MAE_ms",
            ]
        )
        for row in accuracy_rows:
            w.writerow(row)
    print(f"\nSaved: {acc_csv}")

    per_subj_csv = os.path.join(output_dir, "IMU_Synced_per_subject_detection.csv")
    per_subj_df = pd.DataFrame(per_subj_event_rows)
    per_subj_df.to_csv(per_subj_csv, index=False)
    print(f"Saved: {per_subj_csv}")

    alignment_qc_csv = os.path.join(output_dir, "alignment_qc.csv")
    pd.DataFrame(alignment_qc_rows).to_csv(alignment_qc_csv, index=False)
    print(f"Saved: {alignment_qc_csv}")

    # Spatiotemporal summary
    print("\nSPATIOTEMPORAL PARAMETERS BY CONDITION")

    st_rows = []

    for cond_code in cond_list:
        cond_folder = CONDITIONS[cond_code][0]
        print(f"\n  {cond_folder} ({cond_code})")
        print(
            f"  {'Source':<14} {'Stride(s)':>16} {'Step(s)':>16} "
            f"{'Stance(s)':>16} {'Swing(s)':>16} {'Cadence':>14}"
        )

        for src in ["GaitRite", "IMU"]:
            acc = st_accum[cond_code][src]
            stride = np.array(acc["left_stride_time"] + acc["right_stride_time"])
            step = np.array(acc["step_time"])
            stance = np.array(acc["left_stance_time"] + acc["right_stance_time"])
            swing = np.array(acc["left_swing_time"] + acc["right_swing_time"])

            def fmt(arr):
                if len(arr) < 2:
                    return "---"
                return f"{np.mean(arr):.3f}+/-{np.std(arr, ddof=1):.3f}"

            cad_str = "---"
            cad_mean = cad_sd = None
            if len(step) > 1:
                cad = 60.0 / step
                cad_mean = round(float(np.mean(cad)), 1)
                cad_sd = round(float(np.std(cad, ddof=1)), 1)
                cad_str = f"{cad_mean}+/-{cad_sd}"

            print(
                f"  {src:<14} {fmt(stride):>16} {fmt(step):>16} "
                f"{fmt(stance):>16} {fmt(swing):>16} {cad_str:>14}"
            )

            stride_s = summarize_array(stride)
            step_s = summarize_array(step)
            stance_s = summarize_array(stance)
            swing_s = summarize_array(swing)

            st_rows.append(
                [
                    cond_folder,
                    src,
                    stride_s["mean"],
                    stride_s["sd"],
                    stride_s["n"],
                    step_s["mean"],
                    step_s["sd"],
                    step_s["n"],
                    stance_s["mean"],
                    stance_s["sd"],
                    stance_s["n"],
                    swing_s["mean"],
                    swing_s["sd"],
                    swing_s["n"],
                    cad_mean,
                    cad_sd,
                ]
            )

    st_csv = os.path.join(output_dir, "IMU_Synced_spatiotemporal_summary.csv")
    with open(st_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "Condition",
                "Source",
                "Stride_mean",
                "Stride_sd",
                "Stride_n",
                "Step_mean",
                "Step_sd",
                "Step_n",
                "Stance_mean",
                "Stance_sd",
                "Stance_n",
                "Swing_mean",
                "Swing_sd",
                "Swing_n",
                "Cadence_mean",
                "Cadence_sd",
            ]
        )
        for row in st_rows:
            w.writerow([v if v is not None else "" for v in row])
    print(f"Saved: {st_csv}")

    # Stride length summary
    print("\nZUPT STRIDE LENGTH BY CONDITION")

    sl_rows = []

    for cond_code in cond_list:
        cond_folder = CONDITIONS[cond_code][0]
        zupt = np.array(
            sl_accum[cond_code]["left_zupt"] + sl_accum[cond_code]["right_zupt"]
        )
        if len(zupt) > 1:
            print(
                f"  {cond_folder:<14} "
                f"{np.mean(zupt):.3f}+/-{np.std(zupt, ddof=1):.3f} m (n={len(zupt)})"
            )
        else:
            print(f"  {cond_folder:<14} ---")

        zupt_s = summarize_array(zupt)
        sl_rows.append([cond_folder, zupt_s["mean"], zupt_s["sd"], zupt_s["n"]])

    sl_csv = os.path.join(output_dir, "IMU_Synced_stride_length_summary.csv")
    with open(sl_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Condition", "ZUPT_mean", "ZUPT_sd", "ZUPT_n"])
        for row in sl_rows:
            w.writerow([v if v is not None else "" for v in row])
    print(f"Saved: {sl_csv}")

    json_path = os.path.join(output_dir, "IMU_Synced_results.json")
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"Saved: {json_path}")

    if stride_level_rows:
        stride_df = pd.DataFrame(stride_level_rows)
        stride_csv = os.path.join(output_dir, "stride_level_pairs.csv")
        stride_df.to_csv(stride_csv, index=False)
        print(f"\nStride-level pairs: {len(stride_df)} rows")
        errors = stride_df["zupt_sl_m"] - stride_df["gr_sl_m"]
        print(
            f"  Error: bias={errors.mean():.4f} m, "
            f"RMSE={np.sqrt((errors**2).mean()):.4f} m"
        )
        print(f"Saved: {stride_csv}")

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == "__main__":
    main()
