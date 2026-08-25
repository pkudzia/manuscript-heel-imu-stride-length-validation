"""GaitRite parsing, GR Pulse synchronization, and IMU/GaitRite time alignment.

Loads GR Pulse walk boundaries with TDMS/IMU clock drift correction,
splits raw GaitRite exports into OUT/BACK walks, and finds the optimal
time offset that aligns IMU heel strikes with the GaitRite reference.
"""

import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


@dataclass
class Config:
    """Signal-processing and GaitRite parsing settings shared by the loaders below."""

    # IMU settings
    sampling_freq: int = 256
    bandpass_low: float = 0.5
    bandpass_high: float = 10.0
    filter_order: int = 4
    right_gyro_z_col: int = 6
    left_gyro_z_col: int = 34
    invert_left_gyro: bool = True
    invert_right_gyro: bool = False

    # GaitRite column names (note: 'Toe Off ' has a trailing space in the source data)
    gr_col_foot: str = "Left/Right Foot"
    gr_col_heel_on: str = "Heel On"
    gr_col_toe_off: str = "Toe Off "
    gr_left_val: int = 0
    gr_right_val: int = 1

    # Alignment search options
    allow_lr_swap: bool = False


CONFIG = Config()


def calibrate_tdms_to_imu_drift(
    pulse_time: np.ndarray,
    pulse: np.ndarray,
    imu_time: np.ndarray,
    imu_signal: np.ndarray,
    fs: int = None,
) -> Tuple[float, float]:
    """Fit a linear drift model between the TDMS (GR Pulse) and IMU clocks.

    Matches GR Pulse rising edges to the nearest IMU activity onset (via a
    windowed RMS threshold) and fits offset = A + B * tdms_time so that
    downstream alignment can correct for clock drift over long sessions.

    Returns:
        A, B: Drift model coefficients.
    """
    if fs is None:
        fs = CONFIG.sampling_freq

    threshold = (pulse.max() + pulse.min()) / 2
    pulse_high = pulse > threshold
    rising_edges = np.where(np.diff(pulse_high.astype(int)) == 1)[0]

    # RMS envelope for activity detection, computed via convolution
    window = int(0.5 * fs)
    squared_signal = imu_signal**2
    kernel = np.ones(window) / window
    rms = np.sqrt(np.convolve(squared_signal, kernel, mode="same"))

    activity_thresh = 1.5

    tdms_times = []
    imu_times_cal = []

    for i in range(min(7, len(rising_edges))):
        gr_start = pulse_time[rising_edges[i]]
        search_start = max(0, int((gr_start - 5) * fs))
        search_end = min(len(rms), int((gr_start + 2) * fs))

        for j in range(search_start, search_end):
            if rms[j] > activity_thresh:
                imu_activity_start = imu_time[j]
                break
        else:
            continue

        tdms_times.append(gr_start)
        imu_times_cal.append(imu_activity_start)

    if len(tdms_times) < 2:
        print("  Warning: Not enough walks for drift calibration, using defaults")
        return 0.0, 0.0

    tdms_times = np.array(tdms_times)
    imu_times_cal = np.array(imu_times_cal)
    offsets = tdms_times - imu_times_cal

    coeffs = np.polyfit(tdms_times, offsets, 1)
    B = coeffs[0]
    A = coeffs[1]

    return A, B


def load_gr_pulse(filepath: str, imu_data: Dict = None) -> List[Dict]:
    """Load GR Pulse and extract walk boundaries, drift-corrected to the IMU clock.

    Args:
        filepath: Path to GR Pulse CSV file.
        imu_data: IMU data dict for auto-calibration (optional).

    Returns:
        List of walk dictionaries with timing info.

    Raises:
        FileNotFoundError: If the GR Pulse file doesn't exist.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"GR Pulse file not found: {filepath}")

    pulse_df = pd.read_csv(filepath)
    pulse_time = pulse_df["time_s"].values
    pulse = pulse_df["GR_Pulse"].values

    if imu_data is not None:
        A_DRIFT, B_DRIFT = calibrate_tdms_to_imu_drift(
            pulse_time, pulse, imu_data["time"], imu_data["right_filt"]
        )
        print(f"  Drift calibration: offset = {A_DRIFT:.3f} + {B_DRIFT:.4f}*t")

        def tdms_to_imu(t):
            return t - (A_DRIFT + B_DRIFT * t)
    else:

        def tdms_to_imu(t):
            return t

    threshold = (pulse.max() + pulse.min()) / 2
    pulse_high = pulse > threshold
    diff = np.diff(pulse_high.astype(int))

    rising_edges = np.where(diff == 1)[0]
    falling_edges = np.where(diff == -1)[0]

    walks = []
    for i, (r, f) in enumerate(zip(rising_edges, falling_edges)):
        tdms_start = pulse_time[r]
        tdms_end = pulse_time[f]
        walks.append(
            {
                "tdms_walk_num": i + 1,
                "tdms_start": tdms_start,
                "tdms_end": tdms_end,
                "imu_start": tdms_to_imu(tdms_start),
                "imu_end": tdms_to_imu(tdms_end),
                "duration": tdms_end - tdms_start,
            }
        )

    print(f"Loaded GR Pulse: {len(walks)} walks detected")
    return walks


def split_gaitrite_data(filepath: str) -> Dict:
    """Split a raw GaitRite export into OUT and BACK walks.

    Methods (in order of preference):
    1. Natural gap (>10s) in time-sorted data separates OUT and BACK.
    2. Detect where timestamps "jump back" in raw row order (indicates a new walk).
    3. Fall back to splitting at the middle row.

    Raises:
        FileNotFoundError: If the GaitRite file doesn't exist.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"GaitRite file not found: {filepath}")

    df = pd.read_excel(filepath)
    df_clean = df.dropna(
        subset=[CONFIG.gr_col_foot, CONFIG.gr_col_heel_on, CONFIG.gr_col_toe_off]
    )
    df_sorted = df_clean.sort_values(CONFIG.gr_col_heel_on).reset_index(drop=True)

    all_hs = df_sorted[CONFIG.gr_col_heel_on].values
    gaps = np.diff(all_hs)

    if gaps.max() > 10:
        # Natural gap in timestamps (Trial 01 typically)
        split_idx = np.argmax(gaps) + 1
        out_df = df_sorted.iloc[:split_idx]
        back_df = df_sorted.iloc[split_idx:]
        split_method = "natural_gap"
    else:
        # GaitRite exports can interleave walks in raw row order (Walk1
        # events followed by Walk2 events with no time-sorted gap between
        # them). Detect the backward jump that marks the new walk.
        raw_times = df_clean[CONFIG.gr_col_heel_on].values
        split_idx = None

        for i in range(1, len(raw_times)):
            if raw_times[i] < raw_times[i - 1] - 1.0:
                split_idx = i
                break

        if split_idx is not None:
            block1 = df_clean.iloc[:split_idx]
            block2 = df_clean.iloc[split_idx:]
            split_method = "time_jump"
        else:
            mid_idx = len(df_clean) // 2
            block1 = df_clean.iloc[:mid_idx]
            block2 = df_clean.iloc[mid_idx:]
            split_method = "row_order"

        if block1[CONFIG.gr_col_heel_on].min() < block2[CONFIG.gr_col_heel_on].min():
            out_df = block1
            back_df = block2
        else:
            out_df = block2
            back_df = block1

    def extract_events(sub_df):
        left = sub_df[sub_df[CONFIG.gr_col_foot] == CONFIG.gr_left_val]
        right = sub_df[sub_df[CONFIG.gr_col_foot] == CONFIG.gr_right_val]
        return {
            "left_hs": left[CONFIG.gr_col_heel_on].values,
            "right_hs": right[CONFIG.gr_col_heel_on].values,
            "left_to": left[CONFIG.gr_col_toe_off].values,
            "right_to": right[CONFIG.gr_col_toe_off].values,
        }

    return {
        "OUT": extract_events(out_df),
        "BACK": extract_events(back_df),
        "split_method": split_method,
    }


def find_optimal_alignment(
    imu_left_hs: np.ndarray,
    imu_right_hs: np.ndarray,
    gr_left_hs: np.ndarray,
    gr_right_hs: np.ndarray,
    imu_walk_start: float = None,
    search_range: float = 3.0,
    search_step: float = 0.005,
) -> Tuple[float, bool, Dict]:
    """Find the time offset that best aligns IMU heel strikes with GaitRite.

    Searches offsets around an initial estimate and scores each one by the
    mean absolute matching error, normalized by the fraction of GaitRite
    events matched within tolerance, so that spurious low-error/low-match
    offsets don't win over well-matched ones.

    Returns:
        best_offset: Optimal offset to add to GaitRite times.
        swap_lr: True if L/R channels should be swapped.
        alignment_stats: Alignment quality statistics.
    """
    if len(imu_left_hs) == 0 and len(imu_right_hs) == 0:
        return 0.0, False, {"error": "No IMU detections"}

    if len(gr_left_hs) == 0 and len(gr_right_hs) == 0:
        return 0.0, False, {"error": "No GaitRite events"}

    gr_first = min(
        gr_left_hs[0] if len(gr_left_hs) > 0 else np.inf,
        gr_right_hs[0] if len(gr_right_hs) > 0 else np.inf,
    )

    if imu_walk_start is not None:
        initial_offset = imu_walk_start - gr_first
    else:
        imu_first = min(
            imu_left_hs[0] if len(imu_left_hs) > 0 else np.inf,
            imu_right_hs[0] if len(imu_right_hs) > 0 else np.inf,
        )
        initial_offset = imu_first - gr_first

    def compute_alignment_error(offset, swap=False):
        gr_left_aligned = gr_left_hs + offset
        gr_right_aligned = gr_right_hs + offset

        if swap:
            imu_left = imu_right_hs
            imu_right = imu_left_hs
        else:
            imu_left = imu_left_hs
            imu_right = imu_right_hs

        errors = []

        for gr_t in gr_left_aligned:
            if len(imu_left) > 0:
                nearest_idx = np.argmin(np.abs(imu_left - gr_t))
                err = (imu_left[nearest_idx] - gr_t) * 1000
                if abs(err) < 150:
                    errors.append(err)

        for gr_t in gr_right_aligned:
            if len(imu_right) > 0:
                nearest_idx = np.argmin(np.abs(imu_right - gr_t))
                err = (imu_right[nearest_idx] - gr_t) * 1000
                if abs(err) < 150:
                    errors.append(err)

        if len(errors) == 0:
            return np.inf, np.inf, 0, []

        return np.mean(np.abs(errors)), np.std(errors), len(errors), errors

    best_offset = initial_offset
    best_swap = False
    best_score = np.inf
    best_stats = None

    offsets_to_try = np.arange(
        initial_offset - search_range, initial_offset + search_range, search_step
    )

    for offset in offsets_to_try:
        swap_options = [False, True] if CONFIG.allow_lr_swap else [False]
        for swap in swap_options:
            mean_err, std_err, n_matched, errors = compute_alignment_error(offset, swap)

            expected_matches = len(gr_left_hs) + len(gr_right_hs)
            if n_matched < expected_matches * 0.6:
                continue
            match_ratio = n_matched / expected_matches
            score = mean_err / match_ratio

            if score < best_score:
                best_score = score
                best_offset = offset
                best_swap = swap
                best_stats = {
                    "mean_abs_error_ms": mean_err,
                    "std_error_ms": std_err,
                    "n_matched": n_matched,
                    "n_expected": expected_matches,
                    "errors_ms": errors,
                }

    if best_stats is None:
        best_stats = {"error": "No valid alignment found"}
    else:
        best_stats["offset"] = best_offset
        best_stats["swap_lr"] = best_swap
        best_stats["method"] = "pattern_matching"

    return best_offset, best_swap, best_stats
