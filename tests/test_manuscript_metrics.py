"""Independent agreement checks and reproducible public-label resampling."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pingouin as pg
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from _publication_order import publication_order  # noqa: E402
from summarize_manuscript_metrics import icc_a1  # noqa: E402


def test_agreement_formula_matches_independent_library():
    reference = np.array([0.6, 0.8, 0.9, 1.1, 1.3, 1.5])
    estimate = reference + np.array([0.01, -0.03, 0.02, 0.01, -0.02, 0.04])
    long = pd.DataFrame(
        {
            "target": list(range(6)) * 2,
            "rater": ["reference"] * 6 + ["estimate"] * 6,
            "value": np.r_[reference, estimate],
        }
    )
    table = pg.intraclass_corr(
        data=long, targets="target", raters="rater", ratings="value"
    )
    expected = table.loc[table["Type"] == "ICC(A,1)", "ICC"].iloc[0]
    assert icc_a1(reference, estimate) == pytest.approx(expected, abs=1e-12)


def test_resampling_order_is_independent_of_input_row_order(tmp_path):
    order = tmp_path / "order.csv"
    order.write_text("subject\nID_003\nID_001\nID_002\n")
    assert publication_order(["ID_002", "ID_003"], order) == ["ID_003", "ID_002"]
    assert publication_order(["ID_003", "ID_002"], order) == ["ID_003", "ID_002"]


def test_resampling_order_rejects_unknown_or_duplicate_codes(tmp_path):
    order = tmp_path / "order.csv"
    order.write_text("subject\nID_001\n")
    with pytest.raises(ValueError, match="missing"):
        publication_order(["ID_002"], order)
    order.write_text("subject\nID_001\nID_001\n")
    with pytest.raises(ValueError, match="duplicate"):
        publication_order(["ID_001"], order)
