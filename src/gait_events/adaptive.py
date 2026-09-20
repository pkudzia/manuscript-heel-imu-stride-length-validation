"""
Adaptive-threshold gait event detection.

Detects heel strikes and toe-offs using thresholds that adapt to the
signal's amplitude and cadence. No hardcoded amplitude thresholds.

References:
    Adaptation of the sagittal-plane gyroscope negative-peak approach
    (Salarian et al., 2004, IEEE TBME 51(8):1434-1443) with adaptive,
    signal-derived thresholds in the spirit of the adaptive gyroscope
    algorithm of Greene et al. (2010, Med Biol Eng Comput 48(12):1251-1260).
    The autocorrelation stride-frequency estimation, IQR-scaled amplitude
    threshold, and post-peak morphology classification below are this
    project's own combination, not a verbatim reproduction of either.

Algorithm:
    1. Estimate stride frequency via autocorrelation of the gyro Z signal.
    2. Compute signal amplitude statistics (IQR of negative peaks).
    3. Scale peak detection threshold to local amplitude and cadence.
    4. HS: negative peaks followed by stance (near-zero post-peak signal).
    5. TO: negative peaks followed by swing (large positive post-peak signal).
       TO detection is independent of HS (no cascading dependency).

Adaptive thresholds matter here because gait speed, and therefore
gyroscope amplitude, varies widely across older-adult participants and
cognitive-load conditions. Independent TO detection eliminates cascading
error from HS.
"""

from typing import Optional, Tuple

import numpy as np
from scipy.signal import find_peaks


def estimate_stride_frequency(
    sig: np.ndarray, fs: int = 256, freq_range: Tuple[float, float] = (0.8, 2.5)
) -> float:
    """
    Estimate stride frequency via autocorrelation.

    Autocorrelation is more robust to transients and noise than FFT
    for quasi-periodic signals like gait.

    Parameters
    ----------
    sig : np.ndarray
        Gyroscope Z signal [rad/s].
    fs : int
        Sampling frequency [Hz].
    freq_range : tuple
        Expected stride frequency range [Hz] (0.8 Hz = 48 spm,
        2.5 Hz = 150 spm).

    Returns
    -------
    float
        Estimated stride frequency [Hz]. Falls back to 1.5 Hz (90 spm)
        if estimation fails.
    """
    sig_centered = sig - np.mean(sig)

    # Autocorrelation via FFT (efficient for long signals)
    n = len(sig_centered)
    fft_sig = np.fft.fft(sig_centered, n=2 * n)
    acf = np.fft.ifft(fft_sig * np.conj(fft_sig)).real[:n]
    acf /= acf[0]

    lag_min = int(fs / freq_range[1])
    lag_max = min(int(fs / freq_range[0]), n - 1)

    if lag_max <= lag_min or lag_min >= n:
        return 1.5

    acf_window = acf[lag_min : lag_max + 1]
    if len(acf_window) == 0:
        return 1.5

    # First prominent peak in the autocorrelation gives the stride period
    peaks, props = find_peaks(acf_window, height=0.1, distance=int(0.3 * fs))

    if len(peaks) == 0:
        peak_lag = lag_min + np.argmax(acf_window)
    else:
        best_peak = peaks[np.argmax(props["peak_heights"])]
        peak_lag = lag_min + best_peak

    if peak_lag == 0:
        return 1.5

    return fs / peak_lag


def compute_adaptive_threshold(
    sig: np.ndarray,
    stride_freq: float,
    scale_factor: float = 0.35,
    baseline_freq: float = 1.5,
    floor_factor: float = 0.15,
) -> float:
    """
    Compute an amplitude threshold that adapts to signal characteristics.

    Parameters
    ----------
    sig : np.ndarray
        Gyroscope Z signal [rad/s].
    stride_freq : float
        Estimated stride frequency [Hz].
    scale_factor : float
        Threshold as fraction of signal amplitude range.
    baseline_freq : float
        Reference stride frequency for cadence scaling [Hz].
    floor_factor : float
        Minimum threshold as fraction of signal range.

    Returns
    -------
    float
        Adaptive peak detection threshold [rad/s].
    """
    # IQR of negative-peak magnitudes is robust to outliers
    neg_values = sig[sig < 0]
    if len(neg_values) < 10:
        return 1.0

    q25, q75 = np.percentile(np.abs(neg_values), [25, 75])
    iqr = q75 - q25
    threshold = iqr * scale_factor

    # Faster walking produces sharper but sometimes smaller peaks, so
    # scale the threshold down slightly as cadence rises above baseline
    cadence_factor = baseline_freq / max(stride_freq, 0.5)
    cadence_factor = np.clip(cadence_factor, 0.7, 1.5)
    threshold *= cadence_factor

    sig_range = np.max(np.abs(sig))
    floor = sig_range * floor_factor
    threshold = max(threshold, floor)

    # Absolute floor to avoid detecting noise
    threshold = max(threshold, 0.5)

    return threshold


class AdaptiveThresholdDetector:
    """Adaptive threshold gait event detector."""

    def __init__(
        self,
        freq_range: Tuple[float, float] = (0.8, 2.5),
        hs_scale_factor: float = 0.35,
        to_scale_factor: float = 0.30,
        baseline_freq: float = 1.5,
        min_stride: float = 0.4,
        floor_factor: float = 0.15,
        post_peak_window: Tuple[float, float] = (0.08, 0.22),
        to_swing_threshold_factor: float = 0.25,
    ):
        """
        Parameters
        ----------
        freq_range : tuple
            Expected stride frequency range [Hz].
        hs_scale_factor : float
            HS threshold scale relative to signal IQR.
        to_scale_factor : float
            TO threshold scale relative to signal IQR.
        baseline_freq : float
            Reference frequency for cadence normalization [Hz].
        min_stride : float
            Minimum time between consecutive events [s].
        floor_factor : float
            Minimum threshold as fraction of signal range.
        post_peak_window : tuple
            (start, end) of the post-peak classification window [s].
        to_swing_threshold_factor : float
            Fraction of signal range for TO swing validation.
        """
        self.freq_range = freq_range
        self.hs_scale_factor = hs_scale_factor
        self.to_scale_factor = to_scale_factor
        self.baseline_freq = baseline_freq
        self.min_stride = min_stride
        self.floor_factor = floor_factor
        self.post_peak_window = post_peak_window
        self.to_swing_threshold_factor = to_swing_threshold_factor

    def _classify_negative_peaks(
        self, sig: np.ndarray, peak_indices: np.ndarray, fs: int, swing_threshold: float
    ) -> Tuple[list, list]:
        """Classify negative peaks as heel strikes or toe-offs, returning
        (hs_indices, to_indices)."""
        hs_indices = []
        to_indices = []

        for pk_idx in peak_indices:
            post_start = pk_idx + int(self.post_peak_window[0] * fs)
            post_end = min(pk_idx + int(self.post_peak_window[1] * fs), len(sig) - 1)
            if post_end <= post_start:
                continue

            post_seg = sig[post_start:post_end]
            post_mean = np.mean(post_seg)
            post_max = np.max(post_seg)

            # TO: large positive rise (swing phase), using an adaptive
            # swing threshold rather than a fixed value
            if post_max > swing_threshold and post_mean > swing_threshold * 0.5:
                to_indices.append(pk_idx)
            else:
                hs_indices.append(pk_idx)

        return hs_indices, to_indices

    def _enforce_min_spacing(
        self, t: np.ndarray, sig: np.ndarray, indices: list, min_spacing: float
    ) -> list:
        """Keep the event with the most extreme signal value when two
        events fall closer together than min_spacing."""
        if len(indices) == 0:
            return []

        indices = sorted(indices)
        selected = [indices[0]]
        for idx in indices[1:]:
            if idx >= len(t) or selected[-1] >= len(t):
                continue
            if t[idx] - t[selected[-1]] >= min_spacing:
                selected.append(idx)
            else:
                if sig[idx] < sig[selected[-1]]:
                    selected[-1] = idx

        return selected

    def detect_heel_strikes(
        self, t: np.ndarray, sig: np.ndarray, fs: int = 256, **kwargs
    ) -> np.ndarray:
        """Detect heel strikes with adaptive thresholds."""
        stride_freq = estimate_stride_frequency(sig, fs, self.freq_range)

        hs_threshold = compute_adaptive_threshold(
            sig,
            stride_freq,
            self.hs_scale_factor,
            self.baseline_freq,
            self.floor_factor,
        )

        sig_range = np.max(np.abs(sig))
        swing_threshold = sig_range * self.to_swing_threshold_factor
        swing_threshold = max(swing_threshold, 1.0)

        # Minimum peak spacing scales with the estimated stride length
        expected_stride_samples = fs / max(stride_freq, 0.5)
        min_distance = max(int(0.25 * expected_stride_samples), int(0.15 * fs))

        neg_peaks, _ = find_peaks(-sig, height=hs_threshold, distance=min_distance)
        if len(neg_peaks) == 0:
            return np.array([])

        hs_indices, _ = self._classify_negative_peaks(
            sig, neg_peaks, fs, swing_threshold
        )
        if len(hs_indices) == 0:
            return np.array([])

        hs_indices = self._enforce_min_spacing(t, sig, hs_indices, self.min_stride)

        valid = [i for i in hs_indices if i < len(t)]
        return t[np.array(valid)]

    def detect_toe_offs(
        self,
        t: np.ndarray,
        sig: np.ndarray,
        hs_times: Optional[np.ndarray] = None,
        fs: int = 256,
        **kwargs,
    ) -> np.ndarray:
        """
        Detect toe-offs independently using adaptive thresholds.

        Does not depend on hs_times (no cascading dependency). The
        hs_times parameter is accepted for interface compatibility only.
        """
        stride_freq = estimate_stride_frequency(sig, fs, self.freq_range)

        # TO peaks are typically larger and more distinct than HS peaks,
        # so the scale factor is slightly lower
        to_threshold = compute_adaptive_threshold(
            sig,
            stride_freq,
            self.to_scale_factor,
            self.baseline_freq,
            self.floor_factor,
        )

        sig_range = np.max(np.abs(sig))
        swing_threshold = sig_range * self.to_swing_threshold_factor
        swing_threshold = max(swing_threshold, 1.0)

        expected_stride_samples = fs / max(stride_freq, 0.5)
        min_distance = max(int(0.25 * expected_stride_samples), int(0.15 * fs))

        neg_peaks, _ = find_peaks(-sig, height=to_threshold, distance=min_distance)
        if len(neg_peaks) == 0:
            return np.array([])

        _, to_indices = self._classify_negative_peaks(
            sig, neg_peaks, fs, swing_threshold
        )
        if len(to_indices) == 0:
            return np.array([])

        to_indices = self._enforce_min_spacing(t, sig, to_indices, self.min_stride)

        to_times = []
        for pk_idx in to_indices:
            if pk_idx >= len(t):
                continue
            zc_time = self._find_zero_crossing(t, sig, pk_idx, fs)
            to_times.append(zc_time)

        return np.array(sorted(to_times))

    def _find_zero_crossing(
        self, t: np.ndarray, sig: np.ndarray, peak_idx: int, fs: int
    ) -> float:
        """Refine TO timing by finding the zero-crossing after the negative
        peak."""
        max_search = min(peak_idx + int(0.15 * fs), len(sig) - 2)
        for j in range(peak_idx, max_search):
            if sig[j] <= 0 and sig[j + 1] > 0:
                frac = -sig[j] / (sig[j + 1] - sig[j])
                return t[j] + frac * (t[j + 1] - t[j])
        return t[peak_idx]
