#!/usr/bin/env python3
"""Generate manuscript-facing secondary metrics from paired strides."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pingouin as pg
from scipy.stats import mannwhitneyu


def icc_a1(reference, estimate):
    """Two-way random, absolute-agreement, single-measure ICC."""
    values = np.column_stack([reference, estimate]).astype(float)
    n, k = values.shape
    grand = values.mean()
    ss_subject = k * np.square(values.mean(axis=1) - grand).sum()
    ss_method = n * np.square(values.mean(axis=0) - grand).sum()
    residual = values - values.mean(axis=1, keepdims=True) - values.mean(axis=0) + grand
    ss_error = np.square(residual).sum()
    ms_subject = ss_subject / (n - 1)
    ms_method = ss_method / (k - 1)
    ms_error = ss_error / ((n - 1) * (k - 1))
    return (ms_subject - ms_error) / (
        ms_subject + (k - 1) * ms_error + k * (ms_method - ms_error) / n
    )


def agreement_row(data, parameter, condition, reference, estimate):
    data = data[[reference, estimate]].dropna()
    error = data[estimate] - data[reference]
    return {
        "parameter": parameter,
        "condition": condition,
        "gr_mean": data[reference].mean(),
        "gr_sd": data[reference].std(ddof=1),
        "imu_mean": data[estimate].mean(),
        "imu_sd": data[estimate].std(ddof=1),
        "bias": error.mean(),
        "rmse": np.sqrt(np.mean(np.square(error))),
        "icc": icc_a1(data[reference], data[estimate]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", required=True)
    parser.add_argument(
        "--groups",
        default=str(
            Path(__file__).resolve().parents[1] / "metadata" / "participant_groups.csv"
        ),
        help="Anonymized subject,moca,group CSV from Supplementary Table S1",
    )
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    pairs_path = Path(args.pairs)
    output_dir = Path(args.output_dir) if args.output_dir else pairs_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    pairs = pd.read_csv(pairs_path)
    conditions = {"UL", "L2", "L3", "L7"}
    unexpected = set(pairs["condition"].dropna()) - conditions
    if unexpected:
        raise ValueError(
            f"Unexpected conditions in manuscript input: {sorted(unexpected)}"
        )
    if pairs.empty or pairs["condition"].isna().any():
        raise ValueError("Manuscript input must contain labelled stride pairs")
    pairs["gr_speed"] = pairs["gr_sl_m"] / pairs["gr_stride_time_s"]
    pairs["imu_speed"] = pairs["zupt_sl_m"] / pairs["zupt_stride_time_s"]
    pairs["gr_cadence"] = 120 / pairs["gr_stride_time_s"]
    pairs["imu_cadence"] = 120 / pairs["zupt_stride_time_s"]

    rows = []
    parameters = [
        ("Gait speed", "gr_speed", "imu_speed"),
        ("Stride time", "gr_stride_time_s", "zupt_stride_time_s"),
        ("Cadence", "gr_cadence", "imu_cadence"),
    ]
    for parameter, reference, estimate in parameters:
        for condition in ["UL", "L2", "L3", "L7", "Pooled"]:
            data = (
                pairs
                if condition == "Pooled"
                else pairs[pairs["condition"] == condition]
            )
            rows.append(agreement_row(data, parameter, condition, reference, estimate))

    cv_rows = []
    for (subject, condition), data in pairs.groupby(["subject", "condition"]):
        cv_rows.append(
            {
                "subject": subject,
                "condition": condition,
                "gr_cv": 100 * data["gr_sl_m"].std(ddof=1) / data["gr_sl_m"].mean(),
                "imu_cv": 100
                * data["zupt_sl_m"].std(ddof=1)
                / data["zupt_sl_m"].mean(),
            }
        )
    cv = pd.DataFrame(cv_rows)
    for condition in ["UL", "L2", "L3", "L7", "Pooled"]:
        data = cv if condition == "Pooled" else cv[cv["condition"] == condition]
        rows.append(
            agreement_row(data, "Stride-length CV", condition, "gr_cv", "imu_cv")
        )

    metrics = pd.DataFrame(rows)
    metrics.to_csv(output_dir / "manuscript_spatiotemporal_metrics.csv", index=False)

    participant_condition = pairs.groupby(["subject", "condition"], as_index=False)[
        "gr_sl_m"
    ].mean()
    complete = participant_condition.pivot(
        index="subject", columns="condition", values="gr_sl_m"
    ).dropna()
    long = complete.reset_index().melt(
        id_vars="subject", var_name="condition", value_name="gr_sl_m"
    )
    stride_length_aov = pg.rm_anova(
        long,
        dv="gr_sl_m",
        within="condition",
        subject="subject",
        detailed=True,
        correction=True,
    )
    stride_length_aov.insert(0, "n_complete", len(complete))
    stride_length_aov.to_csv(output_dir / "gr_stride_length_rm_anova.csv", index=False)

    demographics = pd.read_csv(args.groups)
    if demographics["subject"].duplicated().any():
        raise ValueError("Group metadata contains duplicate study codes")
    if not set(demographics["group"]).issubset({"CN", "CI"}):
        raise ValueError("Group labels must be CN or CI")
    missing = set(pairs["subject"]) - set(demographics["subject"])
    if missing:
        raise ValueError(f"Missing screening metadata for {sorted(missing)}")

    pairs["abs_error_cm"] = (pairs["zupt_sl_m"] - pairs["gr_sl_m"]).abs() * 100
    participant_mae = pairs.groupby("subject", as_index=False)["abs_error_cm"].mean()
    participant_mae = participant_mae.merge(
        demographics, on="subject", how="left", validate="one_to_one"
    )
    participant_mae.to_csv(output_dir / "participant_mae_by_group.csv", index=False)
    cn = participant_mae.loc[participant_mae["group"] == "CN", "abs_error_cm"]
    mci = participant_mae.loc[participant_mae["group"] == "CI", "abs_error_cm"]
    u_stat, p_value = mannwhitneyu(cn, mci, alternative="two-sided")
    group_summary = (
        participant_mae.groupby("group")["abs_error_cm"]
        .agg(n="size", median="median", minimum="min", maximum="max")
        .reset_index()
    )
    group_summary["mann_whitney_u"] = u_stat
    group_summary["mann_whitney_p"] = p_value
    group_summary.to_csv(output_dir / "participant_mae_group_summary.csv", index=False)

    # Supplementary Table S3 uses the same paired strides, excluding MoCA 17.
    labelled = pairs.merge(demographics, on="subject", validate="many_to_one")
    sensitivity_rows = []
    for rule, data in [
        ("Primary", labelled),
        ("Restricted CI", labelled[labelled["moca"] >= 18]),
    ]:
        for group, cohort in data.groupby("group"):
            errors_cm = 100 * (cohort["zupt_sl_m"] - cohort["gr_sl_m"])
            sensitivity_rows.append(
                {
                    "rule": rule,
                    "group": group,
                    "n_subjects": cohort["subject"].nunique(),
                    "n_pairs": len(cohort),
                    "rmse_cm": np.sqrt(np.mean(errors_cm**2)),
                    "mae_cm": errors_cm.abs().mean(),
                    "bias_cm": errors_cm.mean(),
                }
            )
    pd.DataFrame(sensitivity_rows).to_csv(
        output_dir / "moca_sensitivity_table.csv", index=False
    )

    error = pairs["zupt_sl_m"] - pairs["gr_sl_m"]
    fit = np.polyfit(pairs["gr_sl_m"], pairs["zupt_sl_m"], 1)
    pooled = pd.DataFrame(
        [
            {
                "n_pairs": len(pairs),
                "n_subjects": pairs["subject"].nunique(),
                "mean_gr_stride_length_m": pairs["gr_sl_m"].mean(),
                "within_5_cm_pct": 100 * (error.abs() <= 0.05).mean(),
                "within_10_cm_pct": 100 * (error.abs() <= 0.10).mean(),
                "relative_rmse_pct": 100
                * np.sqrt(np.mean(np.square(error)))
                / pairs["gr_sl_m"].mean(),
                "relative_mae_pct": 100 * error.abs().mean() / pairs["gr_sl_m"].mean(),
                "linear_slope": fit[0],
                "linear_intercept": fit[1],
                "linear_r_squared": np.corrcoef(pairs["gr_sl_m"], pairs["zupt_sl_m"])[
                    0, 1
                ]
                ** 2,
                "participant_mae_min_cm": participant_mae["abs_error_cm"].min(),
                "participant_mae_max_cm": participant_mae["abs_error_cm"].max(),
                "participant_mae_mean_cm": participant_mae["abs_error_cm"].mean(),
                "participant_mae_median_cm": participant_mae["abs_error_cm"].median(),
                "participants_below_4_cm": (participant_mae["abs_error_cm"] < 4).sum(),
            }
        ]
    )
    pooled.to_csv(output_dir / "manuscript_pooled_descriptives.csv", index=False)

    print(metrics.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("\nStride-length condition ANOVA")
    print(stride_length_aov.to_string(index=False))
    print("\nParticipant MAE groups")
    print(group_summary.to_string(index=False))
    print("\nPooled descriptives")
    print(pooled.to_string(index=False))


if __name__ == "__main__":
    main()
