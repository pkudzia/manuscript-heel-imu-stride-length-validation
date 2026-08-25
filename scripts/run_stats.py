#!/usr/bin/env python3
"""Statistical analyses for contralateral ZUPT stride length validation.

Reads stride_level_pairs.csv and produces:
  1. Stride-level error metrics (RMSE, bias, LoA, MAE) by condition
  2. ICC(2,1) between ZUPT and GaitRite per condition
  3. Bland-Altman proportional bias assessment
  4. Mixed-effects model and repeated-measures ANOVA on condition effect
  5. Heteroscedasticity of error vs. gait speed

Writes to the same dated output folder as the input:
  stride_level_error_metrics.csv, icc_results.csv,
  bland_altman_proportional_bias.csv, rm_anova_condition.csv,
  posthoc_condition.csv, statistical_summary.txt
"""

import argparse
import os

import numpy as np
import pandas as pd
import pingouin as pg
import statsmodels.formula.api as smf
from scipy import stats

OUTPUT_DIR = None
STRIDE_CSV = None

COND_ORDER = ["UL", "L2", "L3", "L7"]
COND_LABELS = {
    "UL": "Single-Task",
    "L2": "Serial 2s",
    "L3": "Serial 3s",
    "L7": "Serial 7s",
}


def load_data():
    """Load stride-level pairs for the four reported conditions."""
    df = pd.read_csv(STRIDE_CSV)
    df = df[df["condition"].isin(COND_ORDER)].copy()
    df["error"] = df["zupt_sl_m"] - df["gr_sl_m"]
    df["abs_error"] = df["error"].abs()
    df["mean_sl"] = (df["zupt_sl_m"] + df["gr_sl_m"]) / 2.0
    df["gait_speed_ms"] = df["gr_sl_m"] / df["gr_stride_time_s"]
    return df


def stride_level_metrics(df):
    """Compute error metrics per condition and pooled."""
    results = []

    for cond in COND_ORDER + ["Pooled"]:
        sub = df if cond == "Pooled" else df[df["condition"] == cond]
        n = len(sub)
        if n == 0:
            continue

        err = sub["error"].values
        bias = np.mean(err)
        sd = np.std(err, ddof=1)
        rmse = np.sqrt(np.mean(err**2))
        mae = np.mean(np.abs(err))
        loa_lo = bias - 1.96 * sd
        loa_hi = bias + 1.96 * sd

        subj_means = sub.groupby("subject").agg(
            zupt_mean=("zupt_sl_m", "mean"), gr_mean=("gr_sl_m", "mean")
        )
        subj_err = subj_means["zupt_mean"] - subj_means["gr_mean"]
        subj_rmse = np.sqrt(np.mean(subj_err**2))
        subj_bias = np.mean(subj_err)
        n_subj = len(subj_means)

        label = COND_LABELS.get(cond, cond)
        results.append(
            {
                "condition": cond,
                "label": label,
                "n_strides": n,
                "n_subjects": n_subj,
                "bias_m": round(bias, 4),
                "sd_m": round(sd, 4),
                "rmse_m": round(rmse, 4),
                "mae_m": round(mae, 4),
                "loa_lower_m": round(loa_lo, 4),
                "loa_upper_m": round(loa_hi, 4),
                "subj_level_bias_m": round(subj_bias, 4),
                "subj_level_rmse_m": round(subj_rmse, 4),
            }
        )

    return pd.DataFrame(results)


def compute_icc(df):
    """ICC(2,1) between ZUPT and GaitRite stride length, per condition and pooled."""
    results = []

    for cond in COND_ORDER + ["Pooled"]:
        sub = df if cond == "Pooled" else df[df["condition"] == cond]
        if len(sub) < 10:
            continue

        long = pd.DataFrame(
            {
                "targets": list(range(len(sub))) * 2,
                "raters": ["ZUPT"] * len(sub) + ["GaitRite"] * len(sub),
                "ratings": list(sub["zupt_sl_m"]) + list(sub["gr_sl_m"]),
            }
        )

        icc_res = pg.intraclass_corr(
            data=long, targets="targets", raters="raters", ratings="ratings"
        )
        # ICC(2,1) = two-way random, absolute agreement, single rater.
        # Pingouin labels this 'ICC(A,1)', not 'ICC2'.
        row = icc_res[icc_res["Type"] == "ICC(A,1)"]
        if len(row) == 0:
            continue

        icc_val = row["ICC"].values[0]
        ci_lo = row["CI95"].values[0][0]
        ci_hi = row["CI95"].values[0][1]

        label = COND_LABELS.get(cond, cond)
        results.append(
            {
                "condition": cond,
                "label": label,
                "n_strides": len(sub),
                "icc_2_1": round(icc_val, 4),
                "icc_ci_lower": round(ci_lo, 4),
                "icc_ci_upper": round(ci_hi, 4),
            }
        )

    return pd.DataFrame(results)


def bland_altman_proportional_bias(df):
    """Regress error on the ZUPT/GaitRite mean to test for proportional bias."""
    results = []

    for cond in COND_ORDER + ["Pooled"]:
        sub = df if cond == "Pooled" else df[df["condition"] == cond]
        if len(sub) < 10:
            continue

        slope, intercept, r, p, se = stats.linregress(sub["mean_sl"], sub["error"])
        label = COND_LABELS.get(cond, cond)
        results.append(
            {
                "condition": cond,
                "label": label,
                "slope": round(slope, 4),
                "intercept": round(intercept, 4),
                "r_squared": round(r**2, 4),
                "p_value": round(p, 6),
                "proportional_bias": "Yes" if p < 0.05 else "No",
            }
        )

    return pd.DataFrame(results)


def mixed_effects_condition(df):
    """Mixed-effects model of absolute stride-level error on condition, with
    subject as a random intercept."""
    sub = df[df["condition"].isin(COND_ORDER)].copy()
    sub["condition"] = pd.Categorical(sub["condition"], categories=COND_ORDER)
    model = smf.mixedlm("abs_error ~ C(condition)", sub, groups=sub["subject"])
    return model.fit(reml=True)


def heteroscedasticity_test(df):
    """Regress absolute error on gait speed to check for heteroscedasticity."""
    valid = df.dropna(subset=["gait_speed_ms"])
    if len(valid) < 10:
        return None

    slope, intercept, r, p, se = stats.linregress(
        valid["gait_speed_ms"], valid["abs_error"]
    )
    return {
        "slope": round(slope, 4),
        "intercept": round(intercept, 4),
        "r_squared": round(r**2, 4),
        "p_value": round(p, 6),
        "heteroscedastic": "Yes" if p < 0.05 else "No",
    }


def repeated_measures_anova(df):
    """Repeated-measures ANOVA on subject-level mean absolute error by condition.

    Uses subject-level means (rather than stride-level values) to satisfy
    the independence assumption, restricted to subjects with data in all
    four conditions.
    """
    subj_cond = (
        df.groupby(["subject", "condition"])
        .agg(mean_abs_error=("abs_error", "mean"), n_strides=("abs_error", "count"))
        .reset_index()
    )

    subj_counts = subj_cond.groupby("subject")["condition"].nunique()
    complete_subjects = subj_counts[subj_counts == len(COND_ORDER)].index
    subj_complete = subj_cond[subj_cond["subject"].isin(complete_subjects)]

    if len(complete_subjects) < 5:
        return None, None

    aov = pg.rm_anova(
        data=subj_complete,
        dv="mean_abs_error",
        within="condition",
        subject="subject",
        detailed=True,
    )

    posthoc = pg.pairwise_tests(
        data=subj_complete,
        dv="mean_abs_error",
        within="condition",
        subject="subject",
        padjust="bonf",
    )

    return aov, posthoc


def main():
    global OUTPUT_DIR, STRIDE_CSV
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        required=True,
        help="Path to stride_level_pairs.csv",
    )
    parser.add_argument(
        "--output-dir",
        help="Statistics output directory (default: INPUT/statistics)",
    )
    args = parser.parse_args()

    STRIDE_CSV = os.path.abspath(os.path.expanduser(args.input))
    if not os.path.isfile(STRIDE_CSV):
        parser.error(f"input file not found: {STRIDE_CSV}")
    OUTPUT_DIR = (
        os.path.abspath(os.path.expanduser(args.output_dir))
        if args.output_dir
        else os.path.join(os.path.dirname(STRIDE_CSV), "statistics")
    )

    df = load_data()
    print(
        f"Loaded {len(df)} stride-level pairs from {df['subject'].nunique()} subjects"
    )
    print(f"Conditions: {sorted(df['condition'].unique())}")
    for cond in COND_ORDER:
        n = len(df[df["condition"] == cond])
        ns = df[df["condition"] == cond]["subject"].nunique()
        print(f"  {cond}: {n} strides from {ns} subjects")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    summary_lines = []

    print("\n1. Stride-level error metrics")
    metrics = stride_level_metrics(df)
    print(metrics.to_string(index=False))
    metrics.to_csv(
        os.path.join(OUTPUT_DIR, "stride_level_error_metrics.csv"), index=False
    )
    summary_lines.append("1. STRIDE-LEVEL ERROR METRICS\n")
    summary_lines.append(metrics.to_string(index=False) + "\n\n")

    print("\n2. Intraclass correlation coefficients, ICC(2,1)")
    icc_df = compute_icc(df)
    print(icc_df.to_string(index=False))
    icc_df.to_csv(os.path.join(OUTPUT_DIR, "icc_results.csv"), index=False)
    summary_lines.append("2. INTRACLASS CORRELATION COEFFICIENTS - ICC(2,1)\n")
    summary_lines.append(icc_df.to_string(index=False) + "\n\n")

    print("\n3. Bland-Altman proportional bias test")
    ba_df = bland_altman_proportional_bias(df)
    print(ba_df.to_string(index=False))
    ba_df.to_csv(
        os.path.join(OUTPUT_DIR, "bland_altman_proportional_bias.csv"), index=False
    )
    summary_lines.append("3. BLAND-ALTMAN PROPORTIONAL BIAS TEST\n")
    summary_lines.append(ba_df.to_string(index=False) + "\n\n")

    print("\n4a. Mixed-effects model: condition effect on |error|")
    try:
        me_result = mixed_effects_condition(df)
        print(me_result.summary())
        summary_lines.append("4a. MIXED-EFFECTS MODEL\n")
        summary_lines.append(str(me_result.summary()) + "\n\n")
    except Exception as e:
        print(f"  Mixed-effects model failed: {e}")
        summary_lines.append(f"4a. MIXED-EFFECTS MODEL: Failed ({e})\n\n")

    print("\n4b. Repeated-measures ANOVA (subject-level mean |error|)")
    aov, posthoc = repeated_measures_anova(df)
    if aov is not None:
        print(aov.to_string(index=False))
        print("\nPost-hoc (Bonferroni):")
        print(posthoc.to_string(index=False))
        aov.to_csv(os.path.join(OUTPUT_DIR, "rm_anova_condition.csv"), index=False)
        posthoc.to_csv(os.path.join(OUTPUT_DIR, "posthoc_condition.csv"), index=False)
        summary_lines.append("4b. REPEATED-MEASURES ANOVA\n")
        summary_lines.append(aov.to_string(index=False) + "\n")
        summary_lines.append("Post-hoc (Bonferroni):\n")
        summary_lines.append(posthoc.to_string(index=False) + "\n\n")
    else:
        print("  Insufficient complete-case subjects for RM ANOVA")
        summary_lines.append("4b. RM ANOVA: Insufficient data\n\n")

    print("\n5. Heteroscedasticity: |error| vs. gait speed")
    het = heteroscedasticity_test(df)
    if het:
        for k, v in het.items():
            print(f"  {k}: {v}")
        summary_lines.append("5. HETEROSCEDASTICITY\n")
        for k, v in het.items():
            summary_lines.append(f"  {k}: {v}\n")
        summary_lines.append("\n")

    summary_path = os.path.join(OUTPUT_DIR, "statistical_summary.txt")
    with open(summary_path, "w") as f:
        f.write("STATISTICAL ANALYSIS SUMMARY\n")
        f.write("Contralateral ZUPT Stride Length Validation\n\n")
        f.writelines(summary_lines)
    print(f"\nSaved summary: {summary_path}")


if __name__ == "__main__":
    main()
