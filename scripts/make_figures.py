#!/usr/bin/env python3
"""Generate the two data figures in the manuscript.

Figure 2: stride-length agreement (scatter + Bland-Altman).
Figure 3: per-participant mean absolute error.
Figure 1 is the hand-drawn experimental-setup artwork and is not generated
here. The script reads the primary stride-pair table produced by
run_pipeline.py.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats
from scipy.stats import gaussian_kde

RELEASE_ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = RELEASE_ROOT / "manuscript" / "figures"

COLORS = {
    "secondary": "#b2182b",  # brick red  - mean/bias reference lines
    "grey": "#636363",  # dark grey  - data points and annotations
    "light_grey": "#bdbdbd",  # light grey - zero line
}


def setup_style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "lines.linewidth": 0.8,
            "lines.markersize": 2,
            "axes.linewidth": 0.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "axes.grid": False,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
            "figure.dpi": 150,
            "legend.frameon": True,
            "legend.edgecolor": "#bdbdbd",
            "legend.fancybox": False,
        }
    )


def save_fig(fig, name):
    """Save a figure as PDF (vector, used by LaTeX) and PNG (preview)."""
    fig.savefig(FIG_DIR / f"{name}.pdf", format="pdf")
    fig.savefig(FIG_DIR / f"{name}.png", format="png")
    plt.close(fig)
    print(f"  Saved {name}.pdf / .png")


def figure2_scatter_bland_altman(pairs):
    """Two panels: (a) scatter agreement with a 97.5% KDE density contour,
    (b) Bland-Altman with the full-range proportional-bias regression line.
    Both panels share an enforced square box-aspect."""
    # Authored at 7.0 in wide to match Figure 3, so both scale to
    # \textwidth by the same factor and render at matching text sizes.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.5))

    # Panel (a): scatter with 97.5% KDE contour
    x = pairs["gr_sl_m"].values
    y = pairs["zupt_sl_m"].values

    ax1.scatter(
        x, y, s=4, alpha=0.20, color=COLORS["grey"], edgecolors="none", rasterized=True
    )

    lims = [0.6, 1.8]  # tight to data (strides span 0.65-1.73 m)
    plot_lims = [0.3, 2.1]  # KDE grid stays wide so the contour closes cleanly
    ax1.plot(
        plot_lims,
        plot_lims,
        color=COLORS["grey"],
        linewidth=1.0,
        linestyle=(0, (5, 4)),
        alpha=0.85,
        label="Identity",
        zorder=5,
    )

    xy = np.vstack([x, y])
    kde = gaussian_kde(xy)
    grid_x = np.linspace(plot_lims[0], plot_lims[1], 120)
    grid_y = np.linspace(plot_lims[0], plot_lims[1], 120)
    XX, YY = np.meshgrid(grid_x, grid_y)
    ZZ = kde(np.vstack([XX.ravel(), YY.ravel()])).reshape(XX.shape)

    z_at_pts = kde(xy)
    z_sorted = np.sort(z_at_pts)[::-1]
    cum = np.cumsum(z_sorted) / np.sum(z_sorted)
    level_975 = z_sorted[np.searchsorted(cum, 0.975)]
    # Faint fill inside the 97.5% density region, then the contour outline.
    ax1.contourf(
        XX,
        YY,
        ZZ,
        levels=[level_975, ZZ.max()],
        colors=[COLORS["grey"]],
        alpha=0.12,
        zorder=1,
    )
    ax1.contour(
        XX, YY, ZZ, levels=[level_975], colors=["#000000"], linewidths=[0.9], zorder=10
    )

    slope, intercept, r, p, se = stats.linregress(x, y)

    ax1.set_xlabel("GaitRite stride length (m)", fontsize=12)
    ax1.set_ylabel("IMU stride length (m)", fontsize=12)
    ax1.tick_params(axis="both", labelsize=10)
    ax1.set_xlim(lims)
    ax1.set_ylim(lims)
    ax1.set_xticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
    ax1.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
    ax1.set_box_aspect(1)  # square box so the identity line sits at 45 deg

    try:
        import pingouin as pg

        icc_df = pg.intraclass_corr(
            data=pd.concat(
                [
                    pd.DataFrame(
                        {"targets": range(len(x)), "raters": "GR", "ratings": x}
                    ),
                    pd.DataFrame(
                        {"targets": range(len(x)), "raters": "IMU", "ratings": y}
                    ),
                ],
                ignore_index=True,
            ),
            targets="targets",
            raters="raters",
            ratings="ratings",
        )
        icc_row = icc_df[icc_df["Type"].isin(["ICC2", "ICC(A,1)"])]
        icc_val = icc_row["ICC"].values[0]
    except Exception:
        icc_val = np.nan

    txt = f"n = {len(x):,}\nR² = {r**2:.3f}\nICC = {icc_val:.3f}\nSlope = {slope:.3f}"
    ax1.text(
        0.03,
        0.97,
        txt,
        transform=ax1.transAxes,
        fontsize=10,
        va="top",
        ha="left",
        linespacing=1.35,
    )

    ax1.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color=COLORS["grey"],
                lw=1.0,
                ls=(0, (5, 4)),
                alpha=0.85,
                label="Identity",
            ),
            Line2D([0], [0], color="black", lw=0.9, label="97.5% density region"),
        ],
        loc="lower right",
        fontsize=9,
        frameon=False,
    )
    ax1.set_title("a", fontsize=14, fontweight="bold", loc="left")

    # Panel (b): Bland-Altman
    mean_sl = pairs["mean_sl"].values
    error = pairs["error"].values
    bias = np.mean(error)
    sd = np.std(error, ddof=1)
    loa_lo = bias - 1.96 * sd
    loa_hi = bias + 1.96 * sd

    ax2.scatter(
        mean_sl,
        error,
        s=4,
        alpha=0.20,
        color=COLORS["grey"],
        edgecolors="none",
        rasterized=True,
    )

    # House convention: red solid mean/bias line, black dashed limits.
    ax2.axhline(
        bias,
        color=COLORS["secondary"],
        linewidth=1.2,
        linestyle="-",
        alpha=0.9,
        zorder=5,
    )
    ax2.axhline(
        loa_lo, color="black", linewidth=1.0, linestyle=(0, (5, 4)), alpha=0.9, zorder=5
    )
    ax2.axhline(
        loa_hi, color="black", linewidth=1.0, linestyle=(0, (5, 4)), alpha=0.9, zorder=5
    )
    ax2.axhspan(loa_lo, loa_hi, alpha=0.08, color=COLORS["grey"], zorder=0)

    # Proportional-bias regression line, grey dotted so it reads as a trend
    # distinct from the red mean-bias line.
    prop_slope, prop_int, _, _, _ = stats.linregress(mean_sl, error)
    x_fit2 = np.array(lims)
    ax2.plot(
        x_fit2,
        prop_int + prop_slope * x_fit2,
        color=COLORS["grey"],
        linewidth=1.4,
        linestyle=":",
        alpha=0.95,
        zorder=6,
    )

    ax2.text(
        0.03,
        bias + 0.005,
        f"Bias = {bias:+.3f} m",
        fontsize=10,
        va="bottom",
        ha="left",
        transform=ax2.get_yaxis_transform(),
        color=COLORS["secondary"],
        fontweight="bold",
    )
    ax2.text(
        0.03,
        loa_hi + 0.005,
        f"+1.96 SD = {loa_hi:+.3f} m",
        fontsize=10,
        va="bottom",
        ha="left",
        transform=ax2.get_yaxis_transform(),
        color="black",
        fontweight="bold",
    )
    ax2.text(
        0.03,
        loa_lo - 0.005,
        f"−1.96 SD = {loa_lo:.3f} m",
        fontsize=10,
        va="top",
        ha="left",
        transform=ax2.get_yaxis_transform(),
        color="black",
        fontweight="bold",
    )

    ax2.set_xlabel("Mean stride length (m)", fontsize=12)
    ax2.set_ylabel("IMU − GaitRite (m)", fontsize=12)
    ax2.tick_params(axis="both", labelsize=10)
    ax2.set_xlim(lims)
    ax2.set_ylim(-0.15, 0.15)
    ax2.set_xticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
    ax2.set_yticks([-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15])
    ax2.axhline(0, color=COLORS["light_grey"], linewidth=0.4, zorder=0)
    n_clipped = int(np.sum((error < -0.15) | (error > 0.15)))
    ax2.text(
        0.97,
        0.03,
        f"{n_clipped} strides outside y-axis",
        transform=ax2.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.5,
        color=COLORS["grey"],
    )
    ax2.set_title("b", fontsize=14, fontweight="bold", loc="left")

    for ax in (ax1, ax2):
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)

    ax2.set_box_aspect(1)  # match panel (a) box so the two read as a pair
    fig.tight_layout(w_pad=1.5)
    save_fig(fig, "Figure2")


def figure3_participant_accuracy(pairs):
    """Each dot is one participant's MAE with walk-cluster bootstrap 95%
    confidence intervals, sorted from lowest to highest MAE. Horizontal
    lines show the mean and median participant MAE."""
    pairs = pairs.copy()
    pairs["abs_error"] = np.abs(pairs["error"])
    subj_stats = (
        pairs.groupby("subject")
        .agg(
            mae=("abs_error", "mean"),
            n=("abs_error", "count"),
        )
        .reset_index()
    )
    subj_stats = subj_stats.sort_values("mae").reset_index(drop=True)
    subj_stats["mae_cm"] = subj_stats["mae"] * 100

    # Within-participant walk-cluster bootstrap. Resampling whole walks avoids
    # treating adjacent strides from the same pass as independent observations.
    cluster_cols = ["condition", "trial_num", "walk"]
    missing_cluster_cols = [c for c in cluster_cols if c not in pairs.columns]
    if missing_cluster_cols:
        raise KeyError(f"Missing walk-cluster columns: {missing_cluster_cols}")
    rng = np.random.default_rng(20260806)
    ci_rows = []
    for subject, subject_df in pairs.groupby("subject"):
        clusters = [
            g["abs_error"].to_numpy()
            for _, g in subject_df.groupby(cluster_cols, dropna=False)
        ]
        n_clusters = len(clusters)
        boot_mae = np.empty(2000)
        for b in range(len(boot_mae)):
            draw = rng.integers(0, n_clusters, n_clusters)
            boot_mae[b] = np.mean(np.concatenate([clusters[i] for i in draw]))
        lo, hi = np.percentile(100.0 * boot_mae, [2.5, 97.5])
        ci_rows.append({"subject": subject, "ci_low_cm": lo, "ci_high_cm": hi})
    subj_stats = subj_stats.merge(pd.DataFrame(ci_rows), on="subject", how="left")

    mean_participant_mae = subj_stats["mae_cm"].mean()
    median_participant_mae = subj_stats["mae_cm"].median()

    n_subj = len(subj_stats)
    x = np.arange(n_subj)

    fig, ax = plt.subplots(figsize=(7.0, 3.5))

    # Thin capped 95% CI stems behind the dots
    yerr = np.vstack(
        [
            subj_stats["mae_cm"].values - subj_stats["ci_low_cm"].values,
            subj_stats["ci_high_cm"].values - subj_stats["mae_cm"].values,
        ]
    )
    ax.errorbar(
        x,
        subj_stats["mae_cm"].values,
        yerr=yerr,
        fmt="none",
        ecolor=COLORS["grey"],
        elinewidth=0.9,
        alpha=0.5,
        capsize=2,
        zorder=2,
    )

    ax.scatter(
        x,
        subj_stats["mae_cm"].values,
        s=46,
        color="black",
        edgecolors="white",
        linewidth=0.8,
        zorder=3,
    )

    # House convention: mean = red solid, median = dashed.
    ax.axhline(
        mean_participant_mae,
        color=COLORS["secondary"],
        linewidth=1.0,
        linestyle="-",
        alpha=0.9,
        zorder=4,
    )
    ax.axhline(
        median_participant_mae,
        color="black",
        linewidth=1.0,
        linestyle=(0, (5, 4)),
        alpha=0.9,
        zorder=4,
    )

    ax.set_xlabel("Participant (sorted by MAE)", fontsize=12)
    ax.set_ylabel("Mean absolute error (cm)", fontsize=12)
    ax.tick_params(axis="both", labelsize=10)
    ax.set_xlim(-1, n_subj)
    top_data = subj_stats["ci_high_cm"].max()
    # One tick (0.5 cm) of headroom above the tallest error bar, with a small
    # margin so the top bar never touches a tick line.
    y_top = np.ceil((top_data + 0.2) / 0.5) * 0.5
    ax.set_ylim(0, y_top)
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax.set_xticks([])

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)

    summary = (
        f"n = {n_subj} participants\n"
        f"Total pairs = {subj_stats['n'].sum():,}\n"
        f"Range: {subj_stats['mae_cm'].min():.1f}–"
        f"{subj_stats['mae_cm'].max():.1f} cm"
    )
    ax.text(
        0.02,
        0.97,
        summary,
        transform=ax.transAxes,
        fontsize=10,
        va="top",
        ha="left",
        linespacing=1.35,
    )

    ax.text(
        0.02,
        0.12,
        f"Mean participant MAE = {mean_participant_mae:.1f} cm",
        transform=ax.transAxes,
        fontsize=10,
        color=COLORS["secondary"],
        va="bottom",
        ha="left",
        fontweight="bold",
    )
    ax.text(
        0.02,
        0.03,
        f"Median participant MAE = {median_participant_mae:.1f} cm",
        transform=ax.transAxes,
        fontsize=10,
        color="black",
        va="bottom",
        ha="left",
        fontweight="bold",
    )

    fig.tight_layout()
    save_fig(fig, "Figure3")


def main():
    global FIG_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir", required=True, help="Directory containing stride_level_pairs.csv"
    )
    parser.add_argument(
        "--output-dir", default=str(FIG_DIR), help="Figure output directory"
    )
    args = parser.parse_args()
    data_dir = Path(args.data_dir).resolve()
    FIG_DIR = Path(args.output_dir).resolve()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    setup_style()
    pairs = pd.read_csv(data_dir / "stride_level_pairs.csv")
    pairs = pairs[pairs["condition"].isin(["UL", "L2", "L3", "L7"])].copy()
    pairs["error"] = pairs["zupt_sl_m"] - pairs["gr_sl_m"]
    pairs["mean_sl"] = (pairs["zupt_sl_m"] + pairs["gr_sl_m"]) / 2
    print(f"{len(pairs):,} stride pairs loaded")

    figure2_scatter_bland_altman(pairs)
    figure3_participant_accuracy(pairs)
    print(f"Figures saved to {FIG_DIR}/")


if __name__ == "__main__":
    main()
