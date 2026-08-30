#!/usr/bin/env python3
"""Contralateral vs. unilateral ZUPT stance-detection comparison.

Runs the ZUPT pipeline on all synced participants twice, once with
stance_method='contralateral' (bilateral toe-off validation) and once
with stance_method='threshold' (gyro norm + accel variance only), using
otherwise identical parameters. GaitRite stride pairing uses mutual-
nearest matching within each walk window, matching run_pipeline.py.

Outputs (in the directory selected with ``--output-dir``):
  - contra_vs_uni_stride_pairs.csv   per-stride results for both methods
  - contra_vs_uni_summary.csv        condition-level accuracy comparison
  - contra_vs_uni_stance_stats.csv   stance detection stats
  - contra_vs_uni_report.txt         text summary

Usage:
    python run_contralateral_vs_unilateral.py --base-dir /path/to/Synced-IMU-Data
    python run_contralateral_vs_unilateral.py --conditions UL L2
"""

import argparse
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
    SAMPLING_FREQ,
    imu_column_layout,
    subject_code,
)
from gait_events import AdaptiveThresholdDetector
from stride_length import estimate_stride_length_zupt
from time_sync import find_optimal_alignment, load_gr_pulse, split_gaitrite_data

BASE_DIR = os.environ.get("IMU_SYNCED_BASE_DIR")

# Adaptive detector doubles as the IMU-side anchor for GR alignment.
align_detector = AdaptiveThresholdDetector()


def imu_path(subject, cond_code):
    folder = CONDITIONS[cond_code][0]
    return os.path.join(
        BASE_DIR, "IMU Data", subject, folder, f"{subject}_{cond_code}.csv"
    )


def gr_pulse_path(subject, cond_code):
    for directory in ["Synchronization Data", "Time Sychng Information Data"]:
        base = os.path.join(BASE_DIR, directory, subject)
        for label in [cond_code, CONDITIONS[cond_code][0]]:
            path = os.path.join(base, f"{subject}_{label}_GR_Pulse.csv")
            if os.path.exists(path):
                return path
    return os.path.join(
        BASE_DIR,
        "Synchronization Data",
        subject,
        f"{subject}_{cond_code}_GR_Pulse.csv",
    )


def gaitrite_dir(subject, cond_code):
    folder = CONDITIONS[cond_code][0]
    root = os.path.join(BASE_DIR, "GaitRite Data", subject)
    for candidate in [folder, folder.replace(" ", "")]:
        path = os.path.join(root, candidate)
        if os.path.isdir(path):
            return path
    return os.path.join(root, folder)


def load_imu_csv(filepath):
    """Load IMU CSV, apply gyro-Z polarity fixes, and bandpass-filter."""
    df = pd.read_csv(filepath, skiprows=3, header=None)
    layout = imu_column_layout(df.shape[1])

    time_us = df.iloc[:, 0].values
    time_s = (time_us - time_us[0]) / 1e6

    fs_observed = SAMPLING_FREQ
    if len(time_s) > 1:
        dt = np.median(np.diff(time_s))
        if dt > 0:
            fs_observed = 1.0 / dt

    right_gyro_z = df.iloc[:, layout["right_gyro_z"]].values
    left_gyro_z = df.iloc[:, layout["left_gyro_z"]].values

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

    result["right_accel_3ax"] = df.iloc[:, layout["right_accel"]].values
    result["right_gyro_3ax"] = df.iloc[:, layout["right_gyro"]].values
    result["left_accel_3ax"] = df.iloc[:, layout["left_accel"]].values
    result["left_gyro_3ax"] = df.iloc[:, layout["left_gyro"]].values

    return result


def extract_gaitrite_spatial(filepath):
    """Extract per-walk spatial-temporal data from a GaitRite xlsx.

    Splits OUT/BACK walks where FootFall Object # resets to 1.
    """
    df = pd.read_excel(filepath)
    df = df.dropna(subset=["FootFall Object #"])

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


def load_trial_data(subject, cond_code, imu_data):
    """Load GR Pulse walks and GaitRite trial files for one subject/condition."""
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


def find_subjects(base_dir):
    """Discover subject directories in the synced data folder."""
    subjects = []
    imu_root = os.path.join(base_dir, "IMU Data")
    if not os.path.isdir(imu_root):
        return subjects
    for d in sorted(os.listdir(imu_root)):
        full = os.path.join(imu_root, d)
        if os.path.isdir(full) and not d.startswith("."):
            subjects.append(d)
    return subjects


def run_zupt_both_methods(imu, time_vec, fs):
    """Run ZUPT for both stance methods on both feet.

    Returns a dict keyed by 'contralateral' and 'threshold', each holding
    per-side ZUPT outputs.
    """
    side_config = {
        "left": {"ipsi_filt": imu["left_filt"], "contra_filt": imu["right_filt"]},
        "right": {"ipsi_filt": imu["right_filt"], "contra_filt": imu["left_filt"]},
    }

    results = {}
    for method in ["contralateral", "threshold"]:
        results[method] = {}
        for side in ["left", "right"]:
            accel_key = f"{side}_accel_3ax"
            gyro_key = f"{side}_gyro_3ax"

            if accel_key not in imu or gyro_key not in imu:
                results[method][side] = None
                continue

            try:
                zupt_out = estimate_stride_length_zupt(
                    imu[accel_key],
                    imu[gyro_key],
                    time_vec,
                    fs,
                    ipsi_gyro_z_filt=side_config[side]["ipsi_filt"],
                    contra_gyro_z_filt=side_config[side]["contra_filt"],
                    stance_method=method,
                    bidirectional=True,
                )
                results[method][side] = zupt_out
            except Exception as e:
                print(f"      [{method} {side}] Error: {e}")
                results[method][side] = None

    return results


def filter_to_on_mat(zupt_out, on_mat_windows):
    """Keep only strides whose full HS-to-HS span lies inside a walk window.

    Filtering by midpoint alone would admit strides that straddle the mat
    boundary (e.g. heel-strike before the foot enters the mat), which have
    no valid GaitRite equivalent.
    """
    if zupt_out is None or zupt_out["n_strides"] == 0:
        return {
            "sl": np.array([]),
            "times": np.array([]),
            "st": np.array([]),
        }

    sl = zupt_out["stride_lengths"]
    stride_times_abs = zupt_out["stride_midpoints"]
    st_dur = zupt_out["stride_times"]

    hs_start = stride_times_abs - st_dur / 2.0
    hs_end = stride_times_abs + st_dur / 2.0
    keep = np.zeros(len(sl), dtype=bool)
    for w_start, w_end in on_mat_windows:
        keep |= (hs_start >= w_start) & (hs_end <= w_end)

    mask = keep & np.isfinite(sl)
    return {
        "sl": sl[mask],
        "times": stride_times_abs[mask],
        "st": st_dur[mask],
    }


def icc_absolute_agreement(zupt_vals, gr_vals):
    """Two-way random, absolute-agreement, single-measure ICC(A,1)."""
    ratings = np.column_stack([zupt_vals, gr_vals]).astype(float)
    n, k = ratings.shape
    if n < 2:
        return np.nan
    grand = ratings.mean()
    row_means = ratings.mean(axis=1)
    col_means = ratings.mean(axis=0)
    ms_rows = k * np.sum((row_means - grand) ** 2) / (n - 1)
    ms_cols = n * np.sum((col_means - grand) ** 2) / (k - 1)
    residual = ratings - row_means[:, None] - col_means[None, :] + grand
    ms_error = np.sum(residual**2) / ((n - 1) * (k - 1))
    denominator = ms_rows + (k - 1) * ms_error + k * (ms_cols - ms_error) / n
    return (ms_rows - ms_error) / denominator if denominator != 0 else np.nan


def main():
    global BASE_DIR
    parser = argparse.ArgumentParser(
        description="Contralateral vs. unilateral ZUPT comparison"
    )
    parser.add_argument(
        "--base-dir",
        default=BASE_DIR,
        help="Synced-IMU-Data root (or set IMU_SYNCED_BASE_DIR)",
    )
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--subjects", nargs="+", default=None)
    parser.add_argument("--conditions", nargs="+", default=None)
    args = parser.parse_args()

    if not args.base_dir:
        parser.error("provide --base-dir or set IMU_SYNCED_BASE_DIR")
    BASE_DIR = os.path.abspath(os.path.expanduser(args.base_dir))

    today = date.today().isoformat()
    output_dir = (
        os.path.abspath(os.path.expanduser(args.output_dir))
        if args.output_dir
        else os.path.join(BASE_DIR, "outputs", today, "contra_vs_uni")
    )
    os.makedirs(output_dir, exist_ok=True)

    all_subjects = find_subjects(BASE_DIR)
    subjects = args.subjects if args.subjects else all_subjects
    cond_codes = args.conditions if args.conditions else ["UL", "L2", "L3", "L7"]

    print("Contralateral vs. unilateral ZUPT comparison")
    print(f"Base directory: {BASE_DIR}")
    print(f"Subjects: {len(subjects)}")
    print(f"Conditions: {cond_codes}")
    print(f"Output: {output_dir}\n")

    all_stride_rows = []
    stance_stats_rows = []
    errors_log = []

    for subj in subjects:
        print(f"\nSubject: {subj}")
        subj_label = subject_code(subj)

        for cond_code in cond_codes:
            csv_path = imu_path(subj, cond_code)
            if not os.path.exists(csv_path):
                continue

            cond_folder = CONDITIONS[cond_code][0]
            print(f"  Condition: {cond_code} ({cond_folder})")

            try:
                imu = load_imu_csv(csv_path)
            except Exception as e:
                msg = f"CSV load error for {subj}/{cond_code}: {e}"
                print(f"    {msg}")
                errors_log.append(msg)
                continue

            time_vec = imu["time"]
            fs = imu.get("fs", SAMPLING_FREQ)

            try:
                trials = load_trial_data(subj, cond_code, imu)
            except Exception as e:
                msg = f"Trial load error for {subj}/{cond_code}: {e}"
                print(f"    {msg}")
                errors_log.append(msg)
                continue

            if len(trials) == 0:
                print("    No valid trials found")
                continue

            on_mat_windows = []
            for trial in trials:
                for tdms_walk in [trial["tdms_out"], trial["tdms_back"]]:
                    w_start = tdms_walk.get("imu_start")
                    w_end = tdms_walk.get("imu_end")
                    if w_start is not None and w_end is not None:
                        on_mat_windows.append((w_start, w_end))

            if not on_mat_windows:
                print("    No on-mat windows found")
                continue

            zupt_both = run_zupt_both_methods(imu, time_vec, fs)

            for method in ["contralateral", "threshold"]:
                zupt_with_times = {}
                for side in ["left", "right"]:
                    zout = zupt_both[method].get(side)

                    if zout is not None:
                        stance_stats_rows.append(
                            {
                                "subject": subj_label,
                                "condition": cond_code,
                                "side": side,
                                "stance_method": method,
                                "n_strides": zout["n_strides"],
                                "n_valid": zout.get("n_valid", zout["n_strides"]),
                                "n_stance_phases": len(
                                    zout.get("stride_midpoints", [])
                                ),
                            }
                        )

                    zupt_with_times[side] = filter_to_on_mat(zout, on_mat_windows)

                # Stride-level pairing: time-based matching within each walk
                if (subj_label, cond_code) in EXCLUDE_STRIDE_PAIRS:
                    continue

                for trial in trials:
                    for gr_key, tdms_walk_key in [
                        ("OUT", "tdms_out"),
                        ("BACK", "tdms_back"),
                    ]:
                        if (
                            subj_label,
                            cond_code,
                            trial["trial_num"],
                            gr_key,
                        ) in EXCLUDE_WALKS:
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

                        mask_w = (time_vec >= imu_s - 2.0) & (time_vec <= imu_e + 1.5)
                        t_local = time_vec[mask_w]
                        if len(t_local) < fs:
                            continue
                        try:
                            align_lhs = align_detector.detect_heel_strikes(
                                t_local, imu["left_filt"][mask_w], fs=fs
                            )
                            align_rhs = align_detector.detect_heel_strikes(
                                t_local, imu["right_filt"][mask_w], fs=fs
                            )
                            walk_offset, walk_swap, stats = find_optimal_alignment(
                                align_lhs,
                                align_rhs,
                                gr_raw_hs["left_hs"],
                                gr_raw_hs["right_hs"],
                            )
                            n_matched = int(stats.get("n_matched", 0))
                            n_expected = int(
                                stats.get(
                                    "n_expected",
                                    len(gr_raw_hs["left_hs"])
                                    + len(gr_raw_hs["right_hs"]),
                                )
                            )
                            mean_error = stats.get("mean_abs_error_ms", np.nan)
                            match_ratio = n_matched / n_expected if n_expected else 0.0
                            valid_alignment = (
                                "error" not in stats
                                and np.isfinite(walk_offset)
                                and np.isfinite(mean_error)
                                and n_matched > 0
                                and match_ratio >= 0.60
                            )
                            if not valid_alignment:
                                print(
                                    f"    SYNC QC: excluding {subj} "
                                    f"{cond_code} trial {trial['trial_num']} "
                                    f"{gr_key} from comparison"
                                )
                                continue
                            if walk_swap:
                                print(
                                    f"    WARNING: alignment suggests L/R swap "
                                    f"({subj} {cond_code} trial "
                                    f"{trial['trial_num']} {gr_key}); "
                                    f"check sensor placement."
                                )
                        except (ValueError, IndexError, RuntimeError) as e:
                            print(
                                f"    SYNC QC: excluding {subj} {cond_code} "
                                f"trial {trial['trial_num']} {gr_key} "
                                f"from comparison ({e})"
                            )
                            continue

                        for side in ["left", "right"]:
                            zdata = zupt_with_times[side]
                            if len(zdata["sl"]) == 0:
                                continue

                            in_walk = (zdata["times"] >= imu_s) & (
                                zdata["times"] <= imu_e
                            )
                            z_sl_walk = zdata["sl"][in_walk]
                            z_t_walk = zdata["times"][in_walk]
                            z_st_walk = zdata["st"][in_walk]
                            if len(z_sl_walk) == 0:
                                continue

                            foot_val = 0 if side == "left" else 1
                            foot_mask = gr_foot_arr == foot_val
                            gr_sl_foot = gr_sl_cm[foot_mask] / 100.0
                            gr_st_foot = (
                                gr_st_arr[foot_mask]
                                if len(gr_st_arr) > 0
                                else np.array([])
                            )

                            if len(gr_sl_foot) == 0:
                                continue

                            hs_key = f"{side}_hs"
                            gr_hs_raw = np.sort(np.array(gr_events_raw.get(hs_key, [])))
                            gr_hs = gr_hs_raw + walk_offset
                            n_gr_strides = len(gr_sl_foot)

                            if len(gr_hs) >= 2 and n_gr_strides > 0:
                                n_hs_strides = len(gr_hs) - 1
                                n_use = min(n_gr_strides, n_hs_strides)
                                gr_midpoints = np.array(
                                    [
                                        (gr_hs[i] + gr_hs[i + 1]) / 2.0
                                        for i in range(n_use)
                                    ]
                                )
                                gr_sl_timed = gr_sl_foot[:n_use]
                                gr_st_timed = (
                                    gr_st_foot[:n_use]
                                    if len(gr_st_foot) >= n_use
                                    else np.array([])
                                )
                            else:
                                continue

                            if len(gr_midpoints) == 0:
                                continue

                            # Time bound rejects double-stride artifacts from
                            # missed intermediate heel strikes.
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

                            # No first/last trim: structural completeness is
                            # already guaranteed by ZUPT (bounding stances),
                            # GaitRite (both HS on mat), and the on-mat
                            # filter. Mutual-nearest matching below rejects
                            # any remaining mis-pairs.
                            z_sl_trim = z_sl_walk
                            z_t_trim = z_t_walk

                            # Mutual-nearest matching with cadence-scaled
                            # tolerance (mirrors run_pipeline.py).
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

                                # SL ratio catches stance-detection failures
                                # (double or half-stride errors). Apply the
                                # same stride-time mismatch rule as the
                                # canonical pipeline.
                                gr_sl_val = float(gr_sl_timed[gi])
                                gr_st_val = (
                                    float(gr_st_timed[gi])
                                    if gi < len(gr_st_timed)
                                    else np.nan
                                )
                                zupt_st_val = float(z_st_walk[zi])
                                sl_ratio = (
                                    float(z_sl_trim[zi]) / gr_sl_val
                                    if gr_sl_val > 0
                                    else np.nan
                                )
                                st_ratio = (
                                    zupt_st_val / gr_st_val
                                    if np.isfinite(gr_st_val) and gr_st_val > 0
                                    else np.nan
                                )
                                if not (0.5 <= sl_ratio <= 1.5):
                                    continue
                                if np.isfinite(st_ratio) and not (
                                    0.5 <= st_ratio <= 1.5
                                ):
                                    continue

                                row = {
                                    "subject": subj_label,
                                    "condition": cond_code,
                                    "side": side,
                                    "stance_method": method,
                                    "trial_num": trial["trial_num"],
                                    "walk": gr_key,
                                    "zupt_sl_m": round(float(z_sl_trim[zi]), 4),
                                    "gr_sl_m": round(float(gr_sl_timed[gi]), 4),
                                    "gr_stride_time_s": (
                                        round(float(gr_st_timed[gi]), 4)
                                        if gi < len(gr_st_timed)
                                        else np.nan
                                    ),
                                    "zupt_stride_time_s": round(zupt_st_val, 4),
                                    "zupt_time_s": round(float(z_t_trim[zi]), 4),
                                    "time_diff_s": round(diff_gz, 4),
                                }
                                all_stride_rows.append(row)

                n_method = sum(
                    1
                    for r in all_stride_rows
                    if r["subject"] == subj_label
                    and r["condition"] == cond_code
                    and r["stance_method"] == method
                )
                print(f"    {method}: {n_method} paired strides")

    if not all_stride_rows:
        print("\nNo stride pairs collected. Check data paths and errors above.")
        return

    stride_df = pd.DataFrame(all_stride_rows)
    stride_csv = os.path.join(output_dir, "contra_vs_uni_stride_pairs.csv")
    stride_df.to_csv(stride_csv, index=False)
    print(f"\nSaved: {stride_csv} ({len(stride_df)} rows)")

    # Condition-level summary statistics
    summary_rows = []
    cond_list_with_pooled = ["Pooled"] + cond_codes
    for method in ["contralateral", "threshold"]:
        mdf = stride_df[stride_df["stance_method"] == method]
        for cond in cond_list_with_pooled:
            cdf = (
                mdf[mdf["condition"] != "P"]
                if cond == "Pooled"
                else mdf[mdf["condition"] == cond]
            )
            if len(cdf) == 0:
                continue

            errors = cdf["zupt_sl_m"] - cdf["gr_sl_m"]
            bias = errors.mean()
            sd = errors.std()
            rmse = np.sqrt((errors**2).mean())
            mae = np.abs(errors).mean()
            loa_lower = bias - 1.96 * sd
            loa_upper = bias + 1.96 * sd

            # ICC(A,1): two-way random, absolute agreement, single measure.
            zupt_vals = cdf["zupt_sl_m"].values
            gr_vals = cdf["gr_sl_m"].values
            icc = icc_absolute_agreement(zupt_vals, gr_vals)

            summary_rows.append(
                {
                    "stance_method": method,
                    "condition": cond,
                    "n_subjects": cdf["subject"].nunique(),
                    "n_pairs": len(cdf),
                    "bias_m": round(bias, 4),
                    "sd_m": round(sd, 4),
                    "rmse_m": round(rmse, 4),
                    "mae_m": round(mae, 4),
                    "loa_lower_m": round(loa_lower, 4),
                    "loa_upper_m": round(loa_upper, 4),
                    "icc": round(icc, 3),
                }
            )

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(output_dir, "contra_vs_uni_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"Saved: {summary_csv}")

    print("\nCONTRALATERAL vs. UNILATERAL ZUPT: STRIDE LENGTH ACCURACY")
    print(
        f"{'Method':<16} {'Cond':<8} {'N':>6} {'RMSE':>8} {'Bias':>8} "
        f"{'MAE':>8} {'ICC':>6} {'LoA':>22}"
    )
    for _, row in summary_df.iterrows():
        loa_str = f"[{row['loa_lower_m']:.3f}, {row['loa_upper_m']:.3f}]"
        print(
            f"{row['stance_method']:<16} {row['condition']:<8} "
            f"{row['n_pairs']:>6} {row['rmse_m']:>8.4f} "
            f"{row['bias_m']:>8.4f} {row['mae_m']:>8.4f} "
            f"{row['icc']:>6.3f} {loa_str:>22}"
        )

    contra_pooled = summary_df[
        (summary_df["stance_method"] == "contralateral")
        & (summary_df["condition"] == "Pooled")
    ]
    uni_pooled = summary_df[
        (summary_df["stance_method"] == "threshold")
        & (summary_df["condition"] == "Pooled")
    ]

    if len(contra_pooled) > 0 and len(uni_pooled) > 0:
        c = contra_pooled.iloc[0]
        u = uni_pooled.iloc[0]
        diff = u["rmse_m"] - c["rmse_m"]
        pct = (diff / u["rmse_m"]) * 100 if u["rmse_m"] > 0 else 0

        print("\nKEY RESULT:")
        print(
            f"  Contralateral RMSE: {c['rmse_m']:.4f} m "
            f"(bias={c['bias_m']:.4f}, ICC={c['icc']:.3f}, "
            f"n={c['n_pairs']})"
        )
        print(
            f"  Threshold-only RMSE: {u['rmse_m']:.4f} m "
            f"(bias={u['bias_m']:.4f}, ICC={u['icc']:.3f}, "
            f"n={u['n_pairs']})"
        )
        print(
            f"  Difference: {diff:+.4f} m "
            f"({abs(pct):.1f}% {'improvement' if diff > 0 else 'degradation'})"
        )

    if stance_stats_rows:
        stats_df = pd.DataFrame(stance_stats_rows)
        stats_csv = os.path.join(output_dir, "contra_vs_uni_stance_stats.csv")
        stats_df.to_csv(stats_csv, index=False)
        print(f"\nSaved: {stats_csv}")

        print("\nSTANCE DETECTION STATISTICS")
        for method in ["contralateral", "threshold"]:
            mdf = stats_df[stats_df["stance_method"] == method]
            total = mdf["n_strides"].sum()
            valid = mdf["n_valid"].sum()
            print(f"  {method:<16}: {total:>6} total strides, {valid:>6} valid")

    report_path = os.path.join(output_dir, "contra_vs_uni_report.txt")
    with open(report_path, "w") as f:
        f.write("CONTRALATERAL vs. UNILATERAL ZUPT COMPARISON\n")
        f.write(f"Date: {today}\n")
        f.write(f"Subjects: {len(subjects)}\n")
        f.write(f"Conditions: {cond_codes}\n\n")

        for method in ["contralateral", "threshold"]:
            mdf = stride_df[stride_df["stance_method"] == method]
            if len(mdf) == 0:
                continue
            errors = mdf["zupt_sl_m"] - mdf["gr_sl_m"]
            f.write(f"{method.upper()} (n={len(mdf)} pairs):\n")
            f.write(f"  RMSE = {np.sqrt((errors**2).mean()):.4f} m\n")
            f.write(f"  Bias = {errors.mean():.4f} m\n")
            f.write(f"  MAE  = {np.abs(errors).mean():.4f} m\n")
            f.write(f"  SD   = {errors.std():.4f} m\n\n")

        if errors_log:
            f.write("ERRORS:\n")
            for err in errors_log:
                f.write(f"  - {err}\n")

    print(f"\nSaved: {report_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
