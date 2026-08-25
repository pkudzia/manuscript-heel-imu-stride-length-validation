#!/usr/bin/env python3
"""Validate headline manuscript results against pipeline output tables."""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pingouin as pg

RELEASE_ROOT = Path(__file__).resolve().parent.parent


def assert_close(label, observed, expected, atol=1e-6):
    """Raise a readable assertion when a numeric result has changed."""
    if not np.isclose(observed, expected, atol=atol, rtol=0):
        raise AssertionError(f"{label}: observed {observed!r}, expected {expected!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Directory containing stride_level_pairs.csv and alignment_qc.csv",
    )
    parser.add_argument(
        "--manuscript",
        default=str(RELEASE_ROOT / "manuscript" / "manuscript.tex"),
        help="Manuscript source to check for headline values",
    )
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()
    manifest = json.loads((RELEASE_ROOT / "results_manifest.json").read_text())
    pairs = pd.read_csv(data_dir / "stride_level_pairs.csv")
    alignment = pd.read_csv(data_dir / "alignment_qc.csv")

    pairs = pairs[pairs["condition"].isin(manifest["reported_conditions"])].copy()
    errors = pairs["zupt_sl_m"] - pairs["gr_sl_m"]
    expected_pairs = manifest["stride_pairs"]
    expected_metrics = manifest["pooled_stride_length"]

    assert len(pairs) == expected_pairs["total"]
    assert pairs["subject"].nunique() == manifest["participants"]["analysed"]
    for condition, expected in expected_pairs["by_condition"].items():
        observed = int((pairs["condition"] == condition).sum())
        assert observed == expected, (condition, observed, expected)

    assert_close("pooled RMSE", np.sqrt(np.mean(errors**2)), expected_metrics["rmse_m"])
    assert_close("pooled MAE", np.mean(np.abs(errors)), expected_metrics["mae_m"])
    assert_close("pooled bias", np.mean(errors), expected_metrics["bias_m"])
    long = pd.concat(
        [
            pd.DataFrame(
                {
                    "target": np.arange(len(pairs)),
                    "rater": "GaitRite",
                    "value": pairs["gr_sl_m"].to_numpy(),
                }
            ),
            pd.DataFrame(
                {
                    "target": np.arange(len(pairs)),
                    "rater": "IMU",
                    "value": pairs["zupt_sl_m"].to_numpy(),
                }
            ),
        ],
        ignore_index=True,
    )
    icc_table = pg.intraclass_corr(
        data=long, targets="target", raters="rater", ratings="value"
    )
    pooled_icc = float(icc_table.loc[icc_table["Type"] == "ICC(A,1)", "ICC"].iloc[0])
    assert_close("pooled ICC(2,1)", pooled_icc, expected_metrics["icc_2_1"])

    walk_expected = manifest["walk_alignment"]
    assert len(alignment) == walk_expected["available"]
    valid = alignment["valid"].astype(bool)
    assert int(valid.sum()) == walk_expected["passed"]
    assert int((~valid).sum()) == walk_expected["failed"]
    for condition, expected in walk_expected["failures_by_condition"].items():
        observed = int(((alignment["condition"] == condition) & ~valid).sum())
        assert observed == expected, (condition, observed, expected)

    manuscript = Path(args.manuscript).read_text()
    required = [r"10\{,\}472", "3.9~cm", "0.976", r"1\{,\}669", "95"]
    for token in required:
        if not re.search(token, manuscript):
            raise AssertionError(f"manuscript is missing expected token: {token}")
    stale = [r"8\{,\}884"]
    for token in stale:
        if re.search(token, manuscript):
            raise AssertionError(f"manuscript contains stale result: {token}")

    print(
        "Release validation passed: "
        f"{len(pairs):,} pairs, {valid.sum():,}/{len(alignment):,} aligned walks, "
        f"RMSE {100 * np.sqrt(np.mean(errors**2)):.2f} cm."
    )


if __name__ == "__main__":
    main()
