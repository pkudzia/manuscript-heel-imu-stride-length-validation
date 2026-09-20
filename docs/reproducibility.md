# Reproduction

## Reference analysis

The August 10, 2026 analysis is recorded in [results_manifest.json](../results_manifest.json).
Condition pair counts are 2,544 UL, 2,618 L2, 2,770 L3 and 2,540 L7. There are
1,669 accepted passes and 95 rejected alignments out of 1,764 available passes.

`validate_release.py` recomputes pooled RMSE, MAE, bias and ICC(2,1), verifies
counts and checks manuscript headline values. Numeric tolerance is 0.000001 in
manifest units. This is a whole-cohort regression check, not evidence of validity
on new participants. Manifest `source_outputs` paths record original provenance;
you do not need those local paths to run this release.

`run_stats.py` produces pooled and per-condition agreement, proportional-bias
checks, condition analyses and a text report. `make_figures.py` generates Figures
2 and 3. Figure 1 is supplied experimental-setup artwork.

## Supplementary stance comparison

```bash
python scripts/run_contralateral_vs_unilateral.py \
  --base-dir /path/to/dataset --output-dir outputs/stance_comparison
```

This compares bilateral confirmation with unilateral thresholds. Its rounded
pair table is separate from the primary output.

## Summary and bootstrap outputs

The study period runs from May 2025 to February 2026. After the complete
primary run, create the secondary summary outputs:

```bash
python scripts/summarize_manuscript_metrics.py \
  --pairs outputs/primary/stride_level_pairs.csv --output-dir outputs/summary
python scripts/bootstrap_accuracy.py \
  --pairs outputs/primary/stride_level_pairs.csv \
  --out outputs/summary/bootstrap_ci.csv
```

The summary command writes agreement metrics, repeated-measures stride-length
results, participant MAE group summaries, MoCA sensitivity summaries, and
pooled descriptives. The bootstrap command writes point estimates and 95% CIs
to `bootstrap_ci.csv`. It uses `metadata/participant_groups.csv` through
`scripts/_publication_order.py` to preserve the published subject order. Figure
generation uses the same order file. The ordering correction preserves the
published confidence intervals.

`metadata/participant_groups.csv` contains anonymized study code, screening
score, group labels, and the publication order. It is not a source of raw
recordings or participant linkage.

## Verification

Run `python -m pytest` for the synthetic interface and metric tests. They need
no participant data. Full raw-data reproduction is checked separately, and
participant recordings are never bundled in release assets.
