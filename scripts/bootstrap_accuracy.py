"""
Participant-clustered bootstrap of a primary pipeline stride-pair table.
10,000 iterations. ICC(A,1) computed inline from the two-way random
absolute-agreement single-measures formula (faster than pingouin).
Outputs bootstrap_ci_table4.csv with point estimate + 95% CIs.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from _publication_order import DEFAULT_ORDER, publication_order

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pairs", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--subject-order", default=str(DEFAULT_ORDER))
args = parser.parse_args()
df = pd.read_csv(args.pairs)
df = (
    df[df["condition"].isin(["UL", "L2", "L3", "L7"])]
    .dropna(subset=["zupt_sl_m", "gr_sl_m"])
    .copy()
)
print(f"pairs: {len(df)}, subjects: {df.subject.nunique()}", flush=True)


def icc_A1(g, z):
    """ICC(A,1): two-way random, absolute agreement, single measures."""
    x = np.column_stack([g, z]).astype(float)
    n, k = x.shape
    grand = x.mean()
    sb = ((x.mean(axis=1) - grand) ** 2).sum() * k
    sj = ((x.mean(axis=0) - grand) ** 2).sum() * n
    se = ((x - x.mean(axis=1, keepdims=True) - x.mean(axis=0) + grand) ** 2).sum()
    msb = sb / (n - 1)
    msj = sj / (k - 1)
    mse = se / ((n - 1) * (k - 1))
    return (msb - mse) / (msb + (k - 1) * mse + k * (msj - mse) / n)


def metrics(d):
    err_cm = (d["zupt_sl_m"].values - d["gr_sl_m"].values) * 100.0
    return {
        "rmse_cm": float(np.sqrt(np.mean(err_cm**2))),
        "mae_cm": float(np.mean(np.abs(err_cm))),
        "bias_cm": float(np.mean(err_cm)),
        "icc": float(icc_A1(d["gr_sl_m"].values, d["zupt_sl_m"].values)),
    }


point = metrics(df)
print(f"point: {point}", flush=True)

subjects = publication_order(df.subject.unique(), args.subject_order)
n_subj = len(subjects)
sub_groups = {s: df[df.subject == s].reset_index(drop=True) for s in subjects}
rng = np.random.default_rng(42)

N_ITER = 10000
boot_rmse = np.empty(N_ITER)
boot_mae = np.empty(N_ITER)
boot_bias = np.empty(N_ITER)
boot_icc = np.empty(N_ITER)

for i in range(N_ITER):
    idx = rng.integers(0, n_subj, size=n_subj)
    cohort = pd.concat([sub_groups[subjects[j]] for j in idx], ignore_index=True)
    err_cm = (cohort["zupt_sl_m"].values - cohort["gr_sl_m"].values) * 100.0
    boot_rmse[i] = np.sqrt(np.mean(err_cm**2))
    boot_mae[i] = np.mean(np.abs(err_cm))
    boot_bias[i] = np.mean(err_cm)
    boot_icc[i] = icc_A1(cohort["gr_sl_m"].values, cohort["zupt_sl_m"].values)
    if (i + 1) % 1000 == 0:
        print(f"  iter {i + 1}/{N_ITER}", flush=True)


def ci(arr):
    return float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


rmse_lo, rmse_hi = ci(boot_rmse)
mae_lo, mae_hi = ci(boot_mae)
bias_lo, bias_hi = ci(boot_bias)
icc_lo, icc_hi = ci(boot_icc)

rows = [
    {
        "metric": "rmse_cm",
        "point": point["rmse_cm"],
        "ci_lo": rmse_lo,
        "ci_hi": rmse_hi,
        "n_iter": N_ITER,
    },
    {
        "metric": "mae_cm",
        "point": point["mae_cm"],
        "ci_lo": mae_lo,
        "ci_hi": mae_hi,
        "n_iter": N_ITER,
    },
    {
        "metric": "bias_cm",
        "point": point["bias_cm"],
        "ci_lo": bias_lo,
        "ci_hi": bias_hi,
        "n_iter": N_ITER,
    },
    {
        "metric": "icc",
        "point": point["icc"],
        "ci_lo": icc_lo,
        "ci_hi": icc_hi,
        "n_iter": N_ITER,
    },
]

out = pd.DataFrame(rows)
out_path = Path(args.out)
out_path.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(out_path, index=False)
print(out.to_string(index=False), flush=True)
print(f"saved -> {out_path}", flush=True)
