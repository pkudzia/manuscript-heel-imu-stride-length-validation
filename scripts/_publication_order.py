"""Stable anonymized participant ordering for manuscript bootstrap draws."""

from pathlib import Path

import pandas as pd

DEFAULT_ORDER = Path(__file__).resolve().parents[1] / "metadata/participant_groups.csv"


def publication_order(subjects, order_file=DEFAULT_ORDER):
    """Preserve the published resampling order after de-identification.

    Table S1 preserves the source analysis ordering without exposing source
    identifiers. Sorting the anonymized numeric codes would assign a fixed
    random seed to different participants and change Monte Carlo intervals.
    """
    order = pd.read_csv(order_file)["subject"]
    if order.duplicated().any():
        raise ValueError("Participant-order file contains duplicate study codes")
    selected = set(subjects)
    missing = selected - set(order)
    if missing:
        raise ValueError(f"Participant-order file is missing {sorted(missing)}")
    return [subject for subject in order if subject in selected]
