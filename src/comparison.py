"""Match detected gait events to a GaitRite reference and score the match.

`match_events` pairs detected event times with reference times by
nearest neighbour within a tolerance.
"""

from typing import List, Tuple

import numpy as np


def match_events(
    detected: np.ndarray, reference: np.ndarray, tolerance_s: float
) -> Tuple[List[Tuple[float, float]], List[float], List[float]]:
    """
    Match detected events to reference events via greedy nearest-neighbor.

    Each reference event matches at most one detected event, so a single
    detection can't be double-counted against two reference events.

    Parameters
    ----------
    detected : np.ndarray
        Detected event times [s].
    reference : np.ndarray
        Reference (GaitRite) event times [s].
    tolerance_s : float
        Maximum allowable time difference for a match [s].

    Returns
    -------
    matched : list of (ref_time, det_time) tuples
    missed : list of ref_times with no match (false negatives)
    false_pos : list of det_times with no match (false positives)
    """
    if len(detected) == 0:
        return [], list(reference), []
    if len(reference) == 0:
        return [], [], list(detected)

    used_det = set()
    matched = []
    missed = []

    for ref_t in reference:
        diffs = np.abs(detected - ref_t)
        sorted_idx = np.argsort(diffs)

        found = False
        for idx in sorted_idx:
            if idx in used_det:
                continue
            if diffs[idx] <= tolerance_s:
                matched.append((ref_t, detected[idx]))
                used_det.add(idx)
                found = True
                break
            else:
                break  # remaining candidates are farther away, stop early

        if not found:
            missed.append(ref_t)

    false_pos = [detected[i] for i in range(len(detected)) if i not in used_det]

    return matched, missed, false_pos
