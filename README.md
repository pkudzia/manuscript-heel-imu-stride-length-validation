# Heel-mounted IMU stride-length validation

Private submission-staging repository: <https://github.com/pkudzia/manuscript-heel-imu-stride-length-validation>

Code and manuscript for:

> **Criterion validity of heel-mounted inertial sensors for stride-length estimation during single- and dual-task walking in older adults with and without cognitive impairment**
>
> Pawel Kudzia, Markus von Hacht, Maryam Sheikhi, Drew Commandeur, Marc Klimstra, and Sandra Hundza

The pipeline estimates stride length from two heel-mounted APDM Opal inertial
measurement units and validates the estimates against a GaitRite electronic
walkway. The frozen analysis contains 10,472 matched strides from 32 older
adults across single-task walking and three graded cognitive dual-task
conditions. Pooled stride-length RMSE is 3.86 cm, MAE is 2.50 cm, bias is
+0.26 cm, and ICC(2,1) is 0.976.

## Analysis overview

The method uses adaptive gait-event detection, bilateral stance confirmation,
zero-velocity-aided double integration, adaptive Mahony orientation tracking,
and a retroactive gravity-residual correction. It requires no
participant-specific calibration.

Each GaitRite walk receives a separate alignment decision. A walk enters the
paired validation only when at least 60% of its expected GaitRite heel strikes
match an IMU heel strike within 150 ms. The analysis retained 1,669 of 1,764
available walks. Failed alignment is never replaced with a zero-offset
fallback, and participant-condition records are not removed because their
stride-length error is large.

The frozen headline values are recorded in `results_manifest.json`. The raw
recordings are not included in this release. Participants appear in every
output under de-identified study codes (`ID_NNN`), matching the labels in the
manuscript and supplement.

## Repository layout

```text
src/
  config.py            Sensor layout, shared bounds, documented walk exclusions
  stride_length.py     Bilateral stance detection and ZUPT integration
  time_sync.py         GaitRite parsing, GR Pulse sync, walk alignment
  comparison.py        Event matching
  gait_events/         Adaptive gait-event detector
scripts/
  run_pipeline.py                     Primary analysis and alignment QC
  run_contralateral_vs_unilateral.py  Stance-method comparison
  run_stats.py                        Agreement and condition statistics
  make_figures.py                     Manuscript Figures 2 and 3
  validate_release.py                 Headline-number regression checks
manuscript/                            Audited LaTeX, Word, PDF, bibliography, figures
  manuscript.tex, manuscript.pdf        Current audited LaTeX and PDF files
  manuscript_submission.docx            Current audited Word submission file
  citations/                           Bibliography verification record
  CITATION_AUDIT_2026-08-26.md          Final claim-level citation audit
  SUBMISSION_CHECKLIST.md               Remaining author and upload actions
archive/                               Superseded drafts and working files
submission_upload/                     Files organized for journal upload
results_manifest.json                 Frozen results for the reviewed analysis
```

## Install

Python 3.13 was used for the reviewed analysis.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run the primary pipeline

The data root must contain `IMU Data`, `GaitRite Data`, and
`Time Sychng Information Data` using the study export layout.

```bash
python scripts/run_pipeline.py \
  --base-dir /path/to/Synced-IMU-Data \
  --output-dir /path/to/output/primary
```

The default analysis includes the four reported conditions. Specific
participants or conditions can be selected with `--subjects` and
`--conditions`.

The primary output includes `stride_level_pairs.csv`, `alignment_qc.csv`,
gait-event summaries, and spatiotemporal summaries.

## Statistics, comparison, figures, and validation

```bash
python scripts/run_stats.py \
  --input /path/to/output/primary/stride_level_pairs.csv

python scripts/run_contralateral_vs_unilateral.py \
  --base-dir /path/to/Synced-IMU-Data \
  --output-dir /path/to/output/contra

python scripts/make_figures.py \
  --data-dir /path/to/output/primary \
  --output-dir manuscript/figures

python scripts/validate_release.py \
  --data-dir /path/to/output/primary
```

`validate_release.py` checks participant, walk, and stride counts, pooled
RMSE, MAE, bias, ICC(2,1), per-condition pair counts, and stale manuscript
values against `results_manifest.json`.

## Build the manuscript

```bash
cd manuscript
pdflatex manuscript
bibtex manuscript
pdflatex manuscript
pdflatex manuscript
pdflatex supplementary
pdflatex supplementary
pdflatex cover_letter
```

Build the editable Word files from the repository root:

```bash
python scripts/build_word_doc.py manuscript
python scripts/build_supplementary_word_doc.py manuscript
```

The Word build also requires Pandoc. Built PDF and Word review copies are
included. The main review copies use US Letter pages, line and page numbers,
justified body text, visible paragraph spacing, and intact display-caption
blocks. `Figure1.png` is the experimental-setup artwork and is not generated
by `make_figures.py`.

The LaTeX manuscript uses the official PeerJ `wlpeerj` v1.2 class. The main
PDF retains embedded figures and tables for the LaTeX submission route. The
editable Word copy contains the current Methods figure. PeerJ still receives
the separate artwork files during submission.

## Data and citation

The recordings are APDM Opal IMUs sampled at 256 Hz, with one unit clipped to
the posterolateral heel of each shoe, and a synchronized GaitRite reference.
This private staging repository excludes raw recordings and derived
participant-level output tables. Borealis has reserved DOI
`10.5683/SP4/7YSHYF` for the data deposit. The prospective Data Availability
statement must change to present tense after the dataset and GitHub repository
become public.
