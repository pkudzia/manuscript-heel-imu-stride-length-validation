"""Stride length estimation from a heel-mounted IMU.

The estimator is zero-velocity-aided strapdown integration (ZUPT): detect
stance phases, track foot orientation with an adaptive-gain Mahony filter,
remove gravity, and double-integrate the residual acceleration between
consecutive stance midpoints. Resetting velocity to zero at each stance
boundary bounds the integration drift to a single swing phase.

Entry point: estimate_stride_length_zupt.
"""

from typing import Dict, List, Tuple

import numpy as np

# ZUPT stride length estimation


def _detect_stance_events(
    gyro_z_filt: np.ndarray, time_s: np.ndarray, fs: int = 256
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Detect heel-strike and toe-off times for one leg with the adaptive
    threshold detector.

    Used by the contralateral stance detection method below.

    Args:
        gyro_z_filt: Filtered sagittal gyro-Z signal [rad/s], length N.
        time_s: Time vector [s], length N.
        fs: Sampling frequency [Hz].

    Returns:
        (hs_times, to_times): Arrays of event times in seconds.
    """
    from gait_events import AdaptiveThresholdDetector

    detector = AdaptiveThresholdDetector()
    hs_times = detector.detect_heel_strikes(time_s, gyro_z_filt, fs)
    to_times = detector.detect_toe_offs(time_s, gyro_z_filt, fs=fs)

    return hs_times, to_times


def _detect_stance_contralateral(
    gyro_3ax: np.ndarray,
    accel_3ax: np.ndarray,
    contra_to_times: np.ndarray,
    time_s: np.ndarray,
    fs: int = 256,
    gyro_thresh: float = None,
    accel_thresh: float = None,
    gyro_pct: int = 15,
    accel_pct: int = 25,
) -> np.ndarray:
    """
    Hybrid stance detection: threshold-based candidates, contralateral validation.

    First runs the standard gyro-norm + accel-variance threshold detector to
    find quiet periods (candidate stance windows). Then validates each
    candidate by checking whether a contralateral toe-off event falls within
    or near the window. A contralateral TO confirms ipsilateral single-limb
    support, eliminating false stance triggers from non-gait segments and
    brief signal dips during swing.

    Also recovers missed stance phases: for each contralateral TO not
    already covered by a validated segment, searches for a nearby quiet
    period and adds it.

    Follows the bilateral ZUPT concept of Arens et al. (2021).

    Args:
        gyro_3ax: 3-axis gyroscope [N x 3] in rad/s.
        accel_3ax: 3-axis accelerometer [N x 3] in m/s^2.
        contra_to_times: Contralateral toe-off times [s].
        time_s: Full time vector [s], length N.
        fs: Sampling frequency [Hz].
        gyro_thresh: Gyroscope norm threshold [rad/s]. None = adaptive.
        accel_thresh: Accel variance threshold. None = adaptive.

    Returns:
        Boolean stance mask [N], True = stance.
    """
    N = len(time_s)

    candidate_stance = _detect_stance_zupt(
        gyro_3ax,
        accel_3ax,
        fs,
        gyro_thresh,
        accel_thresh,
        gyro_pct=gyro_pct,
        accel_pct=accel_pct,
    )

    if len(contra_to_times) < 3:
        return candidate_stance

    min_stance_samples = max(int(0.03 * fs), 5)
    merge_gap = max(int(0.075 * fs), 5)
    cand_segs = _find_stance_segments(
        candidate_stance, min_stance_samples, merge_gap_samples=merge_gap
    )

    if len(cand_segs) < 2:
        return candidate_stance

    # A stance segment is valid if a contralateral TO falls within it or
    # within a tolerance window around it.
    tol_s = 0.15  # 150 ms tolerance for event timing misalignment
    tol_samples = int(tol_s * fs)

    validated = np.zeros(N, dtype=bool)
    used_contra_to = set()

    for seg_start, seg_end in cand_segs:
        search_start = time_s[max(0, seg_start - tol_samples)]
        search_end = time_s[min(N - 1, seg_end + tol_samples)]

        matches = np.where(
            (contra_to_times >= search_start) & (contra_to_times <= search_end)
        )[0]

        if len(matches) > 0:
            validated[seg_start : seg_end + 1] = True
            for m in matches:
                used_contra_to.add(m)

    # Recover missed stance phases: for each contralateral TO not already
    # matched, search for a nearby quiet period in the gyro signal and add it.
    gyro_norm = np.linalg.norm(gyro_3ax, axis=1)
    recovery_window_s = 0.3  # search +/- 300 ms around contra TO
    recovery_win = int(recovery_window_s * fs)
    min_quiet_samples = max(int(0.05 * fs), 5)  # 50 ms minimum

    gyro_thresh_recovery = max(np.percentile(gyro_norm, 25), 0.5)

    for i, ct in enumerate(contra_to_times):
        if i in used_contra_to:
            continue

        ct_idx = np.searchsorted(time_s, ct)
        if ct_idx >= N:
            continue

        win_start = max(0, ct_idx - recovery_win)
        win_end = min(N, ct_idx + recovery_win)

        quiet = gyro_norm[win_start:win_end] < gyro_thresh_recovery
        if quiet.sum() < min_quiet_samples:
            continue

        quiet_diff = np.diff(quiet.astype(int))
        q_starts = np.where(quiet_diff == 1)[0] + 1
        q_ends = np.where(quiet_diff == -1)[0] + 1

        if quiet[0]:
            q_starts = np.concatenate([[0], q_starts])
        if quiet[-1]:
            q_ends = np.concatenate([q_ends, [len(quiet)]])

        if len(q_starts) == 0 or len(q_ends) == 0:
            continue

        # Keep the longest contiguous quiet region
        n_pairs = min(len(q_starts), len(q_ends))
        best_len = 0
        best_s, best_e = 0, 0
        for k in range(n_pairs):
            seg_len = q_ends[k] - q_starts[k]
            if seg_len > best_len:
                best_len = seg_len
                best_s = q_starts[k]
                best_e = q_ends[k]

        if best_len >= min_quiet_samples:
            abs_s = win_start + best_s
            abs_e = win_start + best_e
            validated[abs_s:abs_e] = True

    return validated


def _detect_stance_zupt(
    gyro_3ax: np.ndarray,
    accel_3ax: np.ndarray,
    fs: int = 256,
    gyro_thresh: float = None,
    accel_thresh: float = None,
    window_ms: int = 40,
    gyro_pct: int = 15,
    accel_pct: int = 25,
) -> np.ndarray:
    """
    Detect stance phases using a Zero Velocity Update (ZUPT) detector.

    Stance = foot stationary on the ground: gyroscope magnitude near zero
    and accelerometer variance low (only gravity present). Removes
    gyroscope DC bias before thresholding, and uses adaptive percentile
    thresholds when explicit thresholds are not given.

    Args:
        gyro_3ax: 3-axis gyroscope [N x 3] in rad/s.
        accel_3ax: 3-axis accelerometer [N x 3] in m/s^2.
        fs: Sampling frequency [Hz].
        gyro_thresh: Gyroscope norm threshold [rad/s]. None = adaptive.
        accel_thresh: Accelerometer variance threshold [(m/s^2)^2]. None = adaptive.
        window_ms: Sliding window duration [ms].

    Returns:
        Boolean array (True = stance, False = swing) length N.
    """
    N = len(gyro_3ax)
    gyro_norm = np.linalg.norm(gyro_3ax, axis=1)

    win = max(int(fs * window_ms / 1000), 3)
    accel_norm = np.linalg.norm(accel_3ax, axis=1)
    accel_var = np.zeros(N)
    for i in range(N):
        lo = max(0, i - win // 2)
        hi = min(N, i + win // 2 + 1)
        accel_var[i] = np.var(accel_norm[lo:hi])

    if gyro_thresh is None:
        gyro_thresh = max(np.percentile(gyro_norm, gyro_pct), 0.3)
        gyro_thresh = min(gyro_thresh, 2.0)
    if accel_thresh is None:
        accel_thresh = max(np.percentile(accel_var, accel_pct), 0.2)
        accel_thresh = min(accel_thresh, 5.0)

    stance = (gyro_norm < gyro_thresh) & (accel_var < accel_thresh)

    # Remove isolated stance regions shorter than one window
    min_samples = win
    changes = np.diff(stance.astype(int))
    starts = np.where(changes == 1)[0] + 1
    ends = np.where(changes == -1)[0] + 1

    if len(starts) > 0 and len(ends) > 0:
        if ends[0] < starts[0]:
            ends = ends[1:]
        n_pairs = min(len(starts), len(ends))
        for k in range(n_pairs):
            if ends[k] - starts[k] < min_samples:
                stance[starts[k] : ends[k]] = False

    # A "bridge short swing gaps" pass was tried here and reverted 2026-04-28:
    # it caused participant ID_196 to integrate over two real strides as one
    # (~2.4 m strides for a 1.2 m walker). Short-swing merging is instead
    # handled downstream by _find_stance_segments(merge_gap_samples=...).

    # Expect ~40-65% stance during gait; fall back to relaxed percentiles
    # if the detector found almost no stance at all.
    stance_pct = 100.0 * stance.sum() / N
    if stance_pct < 5.0:
        gyro_thresh_r = np.percentile(gyro_norm, 30)
        accel_thresh_r = np.percentile(accel_var, 40)
        stance = (gyro_norm < gyro_thresh_r) & (accel_var < accel_thresh_r)

    return stance


def _quat_multiply(q1, q2):
    """Hamilton quaternion product q1 * q2. Format: [w, x, y, z]."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array(
        [
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        ]
    )


def _quat_to_rotmat(q):
    """Convert a unit quaternion [w, x, y, z] to a 3x3 rotation matrix."""
    w, x, y, z = q / np.linalg.norm(q)
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
        ]
    )


def _quat_from_accel(a):
    """Estimate a tilt quaternion from the gravity vector (no heading information)."""
    a_n = a / np.linalg.norm(a)
    # Gravity in the nav frame points [0, 0, -1]; find the rotation from
    # the measured direction to the nav frame.
    pitch = np.arcsin(np.clip(a_n[0], -1, 1))
    roll = np.arctan2(-a_n[1], -a_n[2])
    cr, sr = np.cos(roll / 2), np.sin(roll / 2)
    cp, sp = np.cos(pitch / 2), np.sin(pitch / 2)
    return np.array([cr * cp, sr * cp, cr * sp, -sr * sp])


def _mahony_update(q, gyro, accel, dt, kp=1.0, ki=0.0, integral_fb=None):
    """
    Mahony AHRS filter update step (Mahony et al., 2008).

    Proportional-integral correction of the gyro using accelerometer-derived
    gravity error.

    Args:
        q: Current quaternion [w, x, y, z].
        gyro: Gyroscope reading [gx, gy, gz] in rad/s (bias-corrected).
        accel: Accelerometer reading [ax, ay, az] in m/s^2.
        dt: Time step (s).
        kp: Proportional gain (higher = faster convergence, more noise).
        ki: Integral gain (compensates residual gyro bias).
        integral_fb: Running integral feedback [3] or None.

    Returns:
        (q_new [4], integral_fb_new [3])
    """
    if integral_fb is None:
        integral_fb = np.zeros(3)

    a_norm = np.linalg.norm(accel)
    if a_norm < 1e-6:
        # No valid accel data; propagate gyro only
        omega_q = np.array([0.0, gyro[0], gyro[1], gyro[2]])
        q_dot = 0.5 * _quat_multiply(q, omega_q)
        q_new = q + q_dot * dt
        return q_new / np.linalg.norm(q_new), integral_fb

    a_hat = accel / a_norm

    # Estimated gravity direction from the current quaternion (third column
    # of the rotation matrix, i.e. rotate [0,0,1] by the conjugate of q)
    w, x, y, z = q
    v_est = np.array(
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            w * w - x * x - y * y + z * z,
        ]
    )

    e = np.cross(a_hat, v_est)  # error between estimated and measured gravity

    integral_fb = integral_fb + ki * e * dt
    gyro_corrected = gyro + kp * e + integral_fb

    omega_q = np.array([0.0, gyro_corrected[0], gyro_corrected[1], gyro_corrected[2]])
    q_dot = 0.5 * _quat_multiply(q, omega_q)
    q_new = q + q_dot * dt
    q_new = q_new / np.linalg.norm(q_new)

    return q_new, integral_fb


def _find_stance_segments(
    stance: np.ndarray, min_dur_samples: int = 10, merge_gap_samples: int = 0
) -> List[Tuple[int, int]]:
    """
    Extract contiguous stance segments from a boolean mask.

    Args:
        stance: Boolean array (True = stance).
        min_dur_samples: Discard segments shorter than this.
        merge_gap_samples: Merge consecutive segments separated by fewer
            than this many samples. A real swing phase lasts 300+ ms; gaps
            shorter than ~75 ms indicate a fragmented stance phase rather
            than a true swing. Set to 0 to disable merging.

    Returns:
        List of (start_idx, end_idx) tuples (inclusive).
    """
    raw_segments = []
    in_st = False
    st_start = 0
    for i in range(len(stance)):
        if stance[i] and not in_st:
            st_start = i
            in_st = True
        elif not stance[i] and in_st:
            if i - st_start >= min_dur_samples:
                raw_segments.append((st_start, i - 1))
            in_st = False
    if in_st and len(stance) - st_start >= min_dur_samples:
        raw_segments.append((st_start, len(stance) - 1))

    if merge_gap_samples <= 0 or len(raw_segments) < 2:
        return raw_segments

    merged = [raw_segments[0]]
    for s, e in raw_segments[1:]:
        prev_s, prev_e = merged[-1]
        gap = s - prev_e - 1
        if gap < merge_gap_samples:
            merged[-1] = (prev_s, e)  # extend previous segment
        else:
            merged.append((s, e))
    return merged


def _estimate_sensor_biases(accel_3ax, gyro_3ax, stance, fs):
    """
    Estimate gyroscope and accelerometer biases from all stance phases.

    Gyro bias: mean gyro during all stance samples (should be zero).
    Accel bias: mean accel during stance minus gravity magnitude along the
    dominant axis (residual after gravity is accounted for).

    Returns:
        (gyro_bias [3], accel_bias [3])
    """
    stance_idx = np.where(stance)[0]
    gyro_bias = np.zeros(3)
    accel_bias = np.zeros(3)

    if len(stance_idx) < 5:
        return gyro_bias, accel_bias

    gyro_bias = np.mean(gyro_3ax[stance_idx], axis=0)

    # Mean stance accel should equal pure gravity; subtract the gravity
    # magnitude from the measured vector to get the residual.
    mean_accel = np.mean(accel_3ax[stance_idx], axis=0)
    g_mag = np.linalg.norm(mean_accel)
    if g_mag > 0:
        g_dir = mean_accel / g_mag
        accel_bias = mean_accel - 9.80665 * g_dir

    return gyro_bias, accel_bias


def _quat_slerp(q0, q1, t):
    """
    Spherical linear interpolation between quaternions q0 and q1.

    Args:
        q0, q1: Unit quaternions [w, x, y, z].
        t: Interpolation parameter (0 = q0, 1 = q1).

    Returns:
        Interpolated unit quaternion [w, x, y, z].
    """
    q0 = q0 / np.linalg.norm(q0)
    q1 = q1 / np.linalg.norm(q1)
    dot = np.dot(q0, q1)
    if dot < 0:  # ensure shortest path
        q1 = -q1
        dot = -dot
    dot = np.clip(dot, -1.0, 1.0)
    if dot > 0.9995:
        # Nearly identical; linear interpolation avoids numerical issues
        result = q0 + t * (q1 - q0)
        return result / np.linalg.norm(result)
    theta = np.arccos(dot)
    sin_theta = np.sin(theta)
    w0 = np.sin((1 - t) * theta) / sin_theta
    w1 = np.sin(t * theta) / sin_theta
    result = w0 * q0 + w1 * q1
    return result / np.linalg.norm(result)


def _zupt_integrate_segment(
    accel_3ax: np.ndarray,
    gyro_3ax: np.ndarray,
    stance: np.ndarray,
    fs: int = 256,
    gyro_bias_override: np.ndarray = None,
    accel_bias_override: np.ndarray = None,
    bidirectional: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Run ZUPT strapdown integration on a data segment.

    Orientation comes from gyro-only quaternion integration with Mahony
    complementary-filter stance corrections, followed by gyro + accel
    bias correction, linear zero-velocity corrections across swing
    phases, and trapezoidal velocity/position integration.

    Args:
        accel_3ax: 3-axis accel [M x 3] in m/s^2.
        gyro_3ax: 3-axis gyro [M x 3] in rad/s.
        stance: Boolean stance mask [M].
        fs: Sampling frequency [Hz].
        gyro_bias_override: Pre-computed gyro bias [3]. When provided,
            skips local bias estimation (used for per-stride reset with
            bout-level bias).
        accel_bias_override: Pre-computed accel bias [3]. Same as above.
        bidirectional: If True, run Mahony forward and backward through
            each swing phase, then SLERP-blend the orientations. Reduces
            mid-swing orientation drift by anchoring to both the preceding
            and following stance phases.

    Returns:
        (velocity [M x 3], position [M x 3], g_nav [3]) in nav frame.
    """
    M = len(accel_3ax)
    dt = 1.0 / fs

    if gyro_bias_override is not None and accel_bias_override is not None:
        gyro_bias = gyro_bias_override
        accel_bias = accel_bias_override
    else:
        gyro_bias, accel_bias = _estimate_sensor_biases(accel_3ax, gyro_3ax, stance, fs)
    accel_corrected = accel_3ax - accel_bias
    gyro_corrected = gyro_3ax - gyro_bias

    # Gravity estimate from the initial stance phase
    stance_idx = np.where(stance)[0]
    if len(stance_idx) >= 3:
        s0 = stance_idx[0]
        s1 = min(s0 + max(int(0.05 * fs), 5), M)
        init_accel = np.mean(accel_corrected[s0:s1], axis=0)
    else:
        init_accel = np.mean(accel_corrected[: max(int(0.05 * fs), 3)], axis=0)

    # Gyro-only orientation: Mahony AHRS (Mahony et al., 2008).
    # Proportional gain kp is higher during stance (gravity observable)
    # and lower during swing (avoids accel artifacts from foot dynamics).
    KP_STANCE = 1.0  # aggressive correction when foot is stationary
    KP_SWING = 0.05  # minimal correction during dynamic motion
    KI = 0.02  # integral gain for residual gyro bias

    q = _quat_from_accel(init_accel)
    q = q / np.linalg.norm(q)
    R0 = _quat_to_rotmat(q)
    g_nav = R0 @ init_accel

    integral_fb = np.zeros(3)
    q_fwd = np.zeros((M, 4))
    q_fwd[0] = q.copy()
    for i in range(1, M):
        kp = KP_STANCE if stance[i] else KP_SWING
        q, integral_fb = _mahony_update(
            q,
            gyro_corrected[i],
            accel_corrected[i],
            dt,
            kp=kp,
            ki=KI,
            integral_fb=integral_fb,
        )
        q_fwd[i] = q.copy()

    accel_nav = np.zeros((M, 3))
    for i in range(M):
        R = _quat_to_rotmat(q_fwd[i])
        accel_nav[i] = R @ accel_corrected[i] - g_nav

    if bidirectional:
        # Retroactive gravity-residual correction.
        #
        # During stance, the foot is stationary, so accel_nav should be
        # zero (gravity fully removed). Any nonzero residual at stance
        # onset reveals tilt drift accumulated during the preceding
        # swing. This residual grows roughly linearly through swing, so
        # the correction is distributed linearly: subtract (t/1) *
        # residual from accel_nav during swing, where t=0 at swing
        # start and t=1 at swing end.
        local_segs = _find_stance_segments(stance, min_dur_samples=3)

        for si in range(len(local_segs) - 1):
            _, end_a = local_segs[si]
            start_b, _ = local_segs[si + 1]
            n_swing = start_b - end_a
            if n_swing <= 1:
                continue

            # Gravity residual at stance B onset. Skip the first ~10 ms
            # (heel-strike impact transient), then average ~20 ms
            # (before Mahony Kp=1.0 has substantially reduced the drift
            # signal).
            skip = max(int(0.010 * fs), 2)  # 10 ms
            win = max(int(0.020 * fs), 3)  # 20 ms
            res_start = min(start_b + skip, M - 1)
            res_end = min(res_start + win, M)
            if res_end <= res_start:
                continue
            g_residual = np.mean(accel_nav[res_start:res_end], axis=0)

            # Linear ramp correction through swing. Scale by 0.5
            # (empirical damping) because drift is not perfectly
            # linear and the residual includes some measurement noise.
            correction_gain = 0.5
            for j in range(1, n_swing + 1):
                idx = end_a + j
                if idx >= M:
                    break
                t = j / n_swing
                accel_nav[idx] -= correction_gain * t * g_residual

    velocity = np.zeros((M, 3))
    for i in range(1, M):
        velocity[i] = velocity[i - 1] + accel_nav[i] * dt

    # Linear ZUPT: distribute accumulated velocity error across each swing
    stance_segs = _find_stance_segments(stance, min_dur_samples=3)
    for si in range(len(stance_segs) - 1):
        _, end_a = stance_segs[si]
        start_b, _ = stance_segs[si + 1]
        n_swing = start_b - end_a
        if n_swing > 0:
            v_s = velocity[end_a].copy()
            v_e = velocity[start_b].copy()
            for j in range(n_swing + 1):
                frac = j / n_swing
                velocity[end_a + j] -= v_s * (1 - frac) + v_e * frac

    for s, e in stance_segs:
        velocity[s : e + 1] = 0.0  # zero velocity during stance

    position = np.zeros((M, 3))
    for i in range(1, M):
        position[i] = position[i - 1] + velocity[i] * dt

    return velocity, position, g_nav


def _split_long_strides(
    stance_mask: np.ndarray,
    stance_segs: List[Tuple[int, int]],
    gyro_3ax: np.ndarray,
    time_s: np.ndarray,
    fs: int = 256,
    max_stride_s: float = 1.8,
) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
    """
    Detect and split strides that are too long (merged double strides).

    Strides longer than max_stride_s are candidates for merged double
    strides caused by a missed stance phase. Searches for a quiet period
    in the gyro norm midway through the stride and inserts a new stance
    segment if found.

    Args:
        stance_mask: Boolean stance mask [N].
        stance_segs: List of (start_idx, end_idx) stance segment tuples.
        gyro_3ax: 3-axis gyroscope [N x 3] in rad/s.
        time_s: Time vector [s], length N.
        fs: Sampling frequency [Hz].
        max_stride_s: Maximum acceptable stride duration [s].

    Returns:
        (updated_stance_mask, updated_stance_segs)
    """
    if len(stance_segs) < 2:
        return stance_mask, stance_segs

    gyro_norm = np.linalg.norm(gyro_3ax, axis=1)
    updated_mask = stance_mask.copy()
    new_segs = list(stance_segs)
    n_splits = 0

    # Use the stance-phase gyro norm values as the quiet-period reference
    stance_gyro_vals = gyro_norm[stance_mask]
    if len(stance_gyro_vals) > 10:
        quiet_thresh = np.percentile(stance_gyro_vals, 85)
    else:
        quiet_thresh = max(np.percentile(gyro_norm, 15), 0.3)
    quiet_thresh = max(quiet_thresh, 0.5)  # floor at 0.5 rad/s

    for i in range(len(stance_segs) - 1):
        mid_a = (stance_segs[i][0] + stance_segs[i][1]) // 2
        mid_b = (stance_segs[i + 1][0] + stance_segs[i + 1][1]) // 2

        stride_time = time_s[mid_b] - time_s[mid_a]

        if stride_time > max_stride_s:
            search_start = stance_segs[i][1] + 1
            search_end = stance_segs[i + 1][0]

            if search_end - search_start < int(0.2 * fs):
                continue

            quiet = gyro_norm[search_start:search_end] < quiet_thresh
            min_quiet = max(int(0.04 * fs), 5)  # 40 ms minimum

            if quiet.sum() < min_quiet:
                relaxed_thresh = np.percentile(gyro_norm[search_start:search_end], 15)
                relaxed_thresh = max(relaxed_thresh, 0.3)
                quiet = gyro_norm[search_start:search_end] < relaxed_thresh
                if quiet.sum() < min_quiet:
                    continue

            quiet_diff = np.diff(quiet.astype(int))
            q_starts = np.where(quiet_diff == 1)[0] + 1
            q_ends = np.where(quiet_diff == -1)[0] + 1

            if quiet[0]:
                q_starts = np.concatenate([[0], q_starts])
            if quiet[-1]:
                q_ends = np.concatenate([q_ends, [len(quiet)]])

            if len(q_starts) == 0 or len(q_ends) == 0:
                continue

            # Keep the longest quiet region
            n_pairs = min(len(q_starts), len(q_ends))
            best_len = 0
            best_s, best_e = 0, 0
            for k in range(n_pairs):
                seg_len = q_ends[k] - q_starts[k]
                if seg_len > best_len:
                    best_len = seg_len
                    best_s = q_starts[k]
                    best_e = q_ends[k]

            if best_len >= min_quiet:
                abs_s = search_start + best_s
                abs_e = search_start + best_e
                new_segs.append((abs_s, abs_e - 1))
                updated_mask[abs_s:abs_e] = True
                n_splits += 1

    if n_splits > 0:
        new_segs = sorted(new_segs, key=lambda x: x[0])

    return updated_mask, new_segs


def estimate_stride_length_zupt(
    accel_3ax: np.ndarray,
    gyro_3ax: np.ndarray,
    time_s: np.ndarray,
    fs: int = 256,
    gyro_thresh: float = None,
    accel_thresh: float = None,
    ipsi_gyro_z_filt: np.ndarray = None,
    contra_gyro_z_filt: np.ndarray = None,
    stance_method: str = "contralateral",
    bidirectional: bool = False,
    gyro_pct: int = 15,
    accel_pct: int = 25,
) -> Dict:
    """
    Self-contained ZUPT stride length estimation from IMU.

    Determines its own gait events from stance-phase detection; no
    external heel-strike input is required. Each stride is defined as the
    interval between the midpoints of consecutive stance phases.

    Processing is done per walking bout (groups of consecutive stance
    phases separated by < 3 s) to limit drift accumulation.

    Algorithm:
        1. Detect stance phases (contralateral or threshold-based).
        2. Group stance phases into walking bouts.
        3. Per bout: gyro-based Mahony orientation tracking, gravity
           removal, double integration with linear ZUPT correction.
        4. Extract horizontal displacement between consecutive stance
           midpoints.

    Args:
        accel_3ax: 3-axis accelerometer [N x 3] in m/s^2.
        gyro_3ax: 3-axis gyroscope [N x 3] in rad/s.
        time_s: Time vector [s], length N.
        fs: Sampling frequency [Hz].
        gyro_thresh: Gyro norm threshold for stance [rad/s]. None = adaptive.
        accel_thresh: Accel variance threshold. None = adaptive.
        ipsi_gyro_z_filt: Filtered sagittal gyro-Z for the ipsilateral leg
            [rad/s]. Required for contralateral stance detection.
        contra_gyro_z_filt: Filtered sagittal gyro-Z for the contralateral
            leg [rad/s]. Required for contralateral stance detection.
        stance_method: 'contralateral' (default) uses bilateral gait
            events following Arens et al. (2021). 'threshold' uses the
            gyro-norm + accel-variance approach alone.

    Returns:
        Dict with:
            'stride_lengths': array of stride lengths [m]
            'stride_times': array of stride durations [s]
            'stride_midpoints': stride midpoint times [s] (HS-equivalent anchor)
            'n_strides': int total strides
            'n_valid': int valid strides (within plausibility bounds)

    References:
        Foxlin, E. (2005). Pedestrian tracking with shoe-mounted
        inertial sensors. IEEE CG&A, 25(6), 38-46.

        Skog, I. et al. (2010). Zero-velocity detection: An algorithm
        evaluation. IEEE TBME, 57(11), 2657-2666.

        Arens, P. et al. (2021). Real-time gait metric estimation for
        everyday gait. Wearable Technologies, 2, e2.
    """
    N = len(accel_3ax)

    use_contralateral = (
        stance_method == "contralateral"
        and ipsi_gyro_z_filt is not None
        and contra_gyro_z_filt is not None
    )

    if use_contralateral:
        _, contra_to = _detect_stance_events(contra_gyro_z_filt, time_s, fs)

        stance = _detect_stance_contralateral(
            gyro_3ax,
            accel_3ax,
            contra_to,
            time_s,
            fs,
            gyro_thresh,
            accel_thresh,
            gyro_pct=gyro_pct,
            accel_pct=accel_pct,
        )

        # Fallback to the threshold detector if the contralateral method
        # was too aggressive. Two triggers:
        #   1. Too few segments overall (< 3)
        #   2. Selectivity < 50% (contralateral rejected more than half of
        #      the threshold-detected stances, which typically means the
        #      contralateral TO detector failed, not that the participant
        #      wasn't walking)
        min_stance_samples = max(int(0.03 * fs), 5)
        merge_gap = max(int(0.075 * fs), 5)
        test_segs = _find_stance_segments(
            stance, min_stance_samples, merge_gap_samples=merge_gap
        )
        thresh_stance = _detect_stance_zupt(
            gyro_3ax,
            accel_3ax,
            fs,
            gyro_thresh,
            accel_thresh,
            gyro_pct=gyro_pct,
            accel_pct=accel_pct,
        )
        thresh_segs = _find_stance_segments(
            thresh_stance, min_stance_samples, merge_gap_samples=merge_gap
        )
        selectivity = len(test_segs) / max(len(thresh_segs), 1)
        if len(test_segs) < 3 or (len(thresh_segs) >= 10 and selectivity < 0.50):
            stance = thresh_stance
    else:
        stance = _detect_stance_zupt(
            gyro_3ax,
            accel_3ax,
            fs,
            gyro_thresh,
            accel_thresh,
            gyro_pct=gyro_pct,
            accel_pct=accel_pct,
        )

    # Find stance segments (each = one foot-flat period). Merge segments
    # separated by < 75 ms to avoid fragmenting a single stance phase into
    # multiple segments: a real swing phase lasts 300+ ms, so anything
    # shorter is noise.
    min_stance_samples = max(int(0.03 * fs), 5)  # at least 30 ms
    merge_gap = max(int(0.075 * fs), 5)  # 75 ms gap threshold
    stance_segs = _find_stance_segments(
        stance, min_stance_samples, merge_gap_samples=merge_gap
    )

    stance, stance_segs = _split_long_strides(
        stance, stance_segs, gyro_3ax, time_s, fs, max_stride_s=1.8
    )

    if len(stance_segs) < 2:
        return {
            "stride_lengths": np.array([]),
            "stride_times": np.array([]),
            "stride_midpoints": np.array([]),
            "n_strides": 0,
            "n_valid": 0,
        }

    # Stance-segment midpoint is used as the ZUPT integration boundary
    # (foot is stationary there, so velocity resets). For the reported
    # stride TIMESTAMP, anchor at stance-segment START (heel-landing
    # equivalent) so the ZUPT stride midpoint matches GaitRite's
    # (HS_i + HS_{i+1})/2 convention -- stance midpoints are offset from
    # HS by ~stance_dur/2 (roughly 300 ms), which eats most of the
    # matching tolerance if used directly.
    stance_mids = np.array([(s + e) // 2 for s, e in stance_segs])
    stance_start_times = np.array([time_s[s] for s, _ in stance_segs])
    if len(stance_start_times) >= 2:
        stride_mid_times = 0.5 * (stance_start_times[:-1] + stance_start_times[1:])
    else:
        stride_mid_times = np.array([])

    # Group into walking bouts (gap > 3 s = new bout)
    max_gap_s = 3.0
    bout_breaks = [0]
    for i in range(1, len(stance_mids)):
        if time_s[stance_segs[i][0]] - time_s[stance_segs[i - 1][1]] > max_gap_s:
            bout_breaks.append(i)
    bout_breaks.append(len(stance_segs))

    all_stride_lengths = []
    all_stride_times = []
    all_mid_times = []

    # Process each bout: estimate biases from all stance phases in the
    # bout, then integrate each stride independently (per-stride reset)
    # so drift does not accumulate across strides.
    for bi in range(len(bout_breaks) - 1):
        bout_start_seg = bout_breaks[bi]
        bout_end_seg = bout_breaks[bi + 1]
        n_segs = bout_end_seg - bout_start_seg

        if n_segs < 2:
            continue

        buf = int(0.2 * fs)
        bout_idx_start = max(0, stance_segs[bout_start_seg][0] - buf)
        bout_idx_end = min(N, stance_segs[bout_end_seg - 1][1] + buf + 1)

        accel_bout = accel_3ax[bout_idx_start:bout_idx_end]
        gyro_bout = gyro_3ax[bout_idx_start:bout_idx_end]
        stance_bout = stance[bout_idx_start:bout_idx_end]

        if len(accel_bout) < 20:
            continue

        gyro_bias, accel_bias = _estimate_sensor_biases(
            accel_bout, gyro_bout, stance_bout, fs
        )

        stance_idx_bout = np.where(stance_bout)[0]
        if len(stance_idx_bout) >= 3:
            s0 = stance_idx_bout[0]
            s1 = min(s0 + max(int(0.05 * fs), 5), len(accel_bout))
            init_accel = np.mean(accel_bout[s0:s1] - accel_bias, axis=0)
        else:
            init_accel = np.mean(
                (accel_bout - accel_bias)[: max(int(0.05 * fs), 3)], axis=0
            )

        g_nav_mag = np.linalg.norm(init_accel)
        g_unit = init_accel / g_nav_mag if g_nav_mag > 0 else np.array([0, 0, 1.0])

        # Per-stride integration: each stride is independent
        for j in range(bout_start_seg, bout_end_seg - 1):
            mid_a = stance_mids[j]
            mid_b = stance_mids[j + 1]

            stride_buf = int(0.05 * fs)
            idx_start = max(0, mid_a - stride_buf)
            idx_end = min(N, mid_b + stride_buf + 1)

            accel_stride = accel_3ax[idx_start:idx_end]
            gyro_stride = gyro_3ax[idx_start:idx_end]
            stance_stride = stance[idx_start:idx_end]

            if len(accel_stride) < 10:
                all_stride_lengths.append(np.nan)
                all_stride_times.append(np.nan)
                all_mid_times.append(stride_mid_times[j])
                continue

            gyro_norm_stride = np.linalg.norm(gyro_stride, axis=1)
            gyro_amplitude = np.ptp(gyro_norm_stride)
            min_gyro_amplitude = 3.0  # rad/s

            if gyro_amplitude < min_gyro_amplitude:
                all_stride_lengths.append(np.nan)
                all_stride_times.append(np.nan)
                all_mid_times.append(stride_mid_times[j])
                continue

            # Integrate this stride independently (position resets to zero)
            vel, pos, g_nav = _zupt_integrate_segment(
                accel_stride,
                gyro_stride,
                stance_stride,
                fs,
                gyro_bias_override=gyro_bias,
                accel_bias_override=accel_bias,
                bidirectional=bidirectional,
            )

            local_mid_a = mid_a - idx_start
            local_mid_b = mid_b - idx_start

            if local_mid_a < 0 or local_mid_b >= len(pos):
                all_stride_lengths.append(np.nan)
                all_stride_times.append(np.nan)
                all_mid_times.append(stride_mid_times[j])
                continue

            disp = pos[local_mid_b] - pos[local_mid_a]
            st = time_s[stance_mids[j + 1]] - time_s[stance_mids[j]]

            # Project onto the horizontal plane (perpendicular to gravity)
            g_unit_local = (
                g_nav / np.linalg.norm(g_nav) if np.linalg.norm(g_nav) > 0 else g_unit
            )
            vert = np.dot(disp, g_unit_local) * g_unit_local
            horiz = disp - vert
            sl = np.linalg.norm(horiz)

            # Validity: 0.3-2.5 m, 0.4-2.0 s. The lower length floor
            # matches GR_SL_MIN_M=0.30 so shuffling/MCI strides are not
            # asymmetrically dropped on the IMU side; the upper time
            # bound admits very slow walkers.
            if np.isfinite(sl) and 0.3 <= sl <= 2.5 and 0.4 <= st <= 2.0:
                all_stride_lengths.append(sl)
            else:
                all_stride_lengths.append(np.nan)
            all_stride_times.append(st)
            all_mid_times.append(stride_mid_times[j])

    stride_lengths = np.array(all_stride_lengths)
    stride_times = np.array(all_stride_times)
    mid_times = np.array(all_mid_times)

    return {
        "stride_lengths": stride_lengths,
        "stride_times": stride_times,
        "stride_midpoints": mid_times,
        "n_strides": len(stride_lengths),
        "n_valid": int(np.sum(np.isfinite(stride_lengths))),
    }
