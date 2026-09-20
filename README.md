# Heel-mounted IMU stride-length validation

Analysis code for the bilateral heel-IMU stride-length validation study.

Pawel Kudzia, Markus von Hacht, Maryam Sheikhi, Drew Commandeur, Marc Klimstra,
and Sandra Hundza.

The pipeline estimates stride length from bilateral heel-mounted APDM Opal IMUs
and compares it with a synchronized GaitRite walkway. The reference analysis
contains **10,472 stride pairs from 32 participants**, with RMSE **3.86 cm**,
MAE **2.50 cm**, bias **+0.26 cm**, and ICC(2,1) **0.976**. Exact regression
targets are in [results_manifest.json](results_manifest.json).

This archive holds the analysis code only. The participant recordings are
published separately through Borealis, and the manuscript is not included.

## Start here

- [Data download](https://doi.org/10.5683/SP4/7YSHYF), Borealis version 1.0.
- [Algorithm and code map](docs/algorithm.md), [input format](docs/data-format.md).
- [Reproduction steps](docs/reproducibility.md), [data access](DATA_ACCESS.md).

## Install

Use **Python 3.13**. Dependencies are pinned to the reproduction environment.
Run these commands in the archive root:

```bash
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate with `.venv\Scripts\activate` instead.

## Prepare the data

Download and extract [Borealis version 1.0](https://doi.org/10.5683/SP4/7YSHYF).
Keep the data outside this tree. Preserve the folders and their contents:

```text
dataset/
  IMU Data/
  GaitRite Data/
  Synchronization Data/
  SHA256SUMS_EXTRACTED_DATA.txt
```

The deposit provides 128 IMU CSV files, 126 synchronization CSV files and 897
GaitRite XLSX files. Verify all 1,151 raw files using the deposited manifest.
On macOS, run this inside the data directory:

```bash
shasum -a 256 -c SHA256SUMS_EXTRACTED_DATA.txt
```

On Linux, use `sha256sum -c SHA256SUMS_EXTRACTED_DATA.txt`.
See [data access](DATA_ACCESS.md) for the two known missing synchronization files.

## Reproduce the primary results

Replace `/path/to/dataset` with the directory containing the three data folders:

```bash
python scripts/run_pipeline.py \
  --base-dir /path/to/dataset --output-dir outputs/primary

python scripts/validate_release.py --data-dir outputs/primary

python scripts/run_stats.py \
  --input outputs/primary/stride_level_pairs.csv

python scripts/make_figures.py \
  --data-dir outputs/primary --output-dir outputs/figures
```

Validation must report **10,472 pairs and 1,669/1,764 aligned walks**.
Run it after a complete cohort analysis. A deliberate subset is useful for
exploration, but is not expected to pass full-cohort reference checks:

```bash
python scripts/run_pipeline.py --base-dir /path/to/dataset \
  --subjects ID_101 --conditions UL --output-dir outputs/example
```

Conditions are `UL`, `L2`, `L3` and `L7`. Invalid selections are rejected.
Use `--help` on each analysis script for its arguments.

## Summary tables

Generate the secondary tables from a completed primary pair table:

```bash
python scripts/summarize_manuscript_metrics.py \
  --pairs outputs/primary/stride_level_pairs.csv --output-dir outputs/summary
python scripts/bootstrap_accuracy.py \
  --pairs outputs/primary/stride_level_pairs.csv \
  --out outputs/summary/bootstrap_ci.csv
```

Both commands read `metadata/participant_groups.csv` through
`scripts/_publication_order.py` to keep the published subject order, which the
fixed-seed bootstrap depends on.

## How it works

1. Read the two heel signals and detect heel strikes and toe-offs.
2. Confirm stance bilaterally and integrate acceleration between stance anchors.
3. Correct clock drift and align each walkway pass.
4. Match eligible same-foot strides within accepted passes.
5. Export paired measurements, alignment decisions and agreement statistics.

The method uses no participant-specific calibration. A pass requires at least
60% of expected walkway heel strikes to match within 150 ms. Failed alignment
is never replaced with a zero-offset fallback. See the
[algorithm guide](docs/algorithm.md).

## Contents

| Location | Purpose |
| --- | --- |
| `src/` | Signal processing, stride integration, synchronization and matching |
| `scripts/` | Analysis commands, verification and figures |
| `tests/` | Synthetic tests that need no participant recordings |
| `docs/` | Algorithm, data format and reproduction instructions |
| `metadata/` | Anonymized group labels and the publication subject order |
| `results_manifest.json` | Frozen primary counts and agreement targets |
| `CITATION.cff` | Versioned software citation |

Generated outputs belong in the ignored `outputs/` directory.

## Development checks

```bash
python -m pip install -r requirements-dev.txt
python -m ruff check src scripts tests
python -m pytest
```

## Citation and licence

Cite the code with [CITATION.cff](CITATION.cff), and cite the
[Borealis data DOI](https://doi.org/10.5683/SP4/7YSHYF) separately.

The code is [MIT licensed](LICENSE). That licence does not relicense the data.
