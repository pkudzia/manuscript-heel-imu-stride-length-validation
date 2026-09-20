# Input data and outputs

Use the extracted [Borealis version 1.0 data](https://doi.org/10.5683/SP4/7YSHYF).
Preserve directories, file names and column labels.

```text
dataset/
  IMU Data/ID_101/Unloaded/ID_101_UL.csv
  GaitRite Data/ID_101/Unloaded/ID_101_Unloaded_Trial_01.xlsx
  Synchronization Data/ID_101/ID_101_UL_GR_Pulse.csv
```

| Code | Folder | Condition |
| --- | --- | --- |
| UL | Unloaded | Single task |
| L2 | Loaded 2 | Serial subtraction by 2 |
| L3 | Loaded 3 | Serial subtraction by 3 |
| L7 | Loaded 7 | Serial subtraction by 7 |

GaitRite folders may also use `Loaded2` without a space. The legacy synchronization
directory name `Time Sychng Information Data` is accepted. Outputs use `ID_NNN` codes.

## IMU CSV

Skip three preamble rows, then read numeric rows without a header. Time is in
microseconds, acceleration in m/s² and angular velocity in rad/s. Nominal rate
is 256 Hz; actual rate is estimated from timestamp differences.

These are **zero-based** column indices. Other APDM fields are not consumed by
the stride estimator.

| Field | Public 29-column export | Original 43-column export |
| --- | --- | --- |
| Time | 0 | 0 |
| Right acceleration x/y/z | 1, 2, 3 | 1, 2, 3 |
| Right gyroscope x/y/z | 4, 5, 6 | 4, 5, 6 |
| Left acceleration x/y/z | 15, 16, 17 | 29, 30, 31 |
| Left gyroscope x/y/z | 18, 19, 20 | 32, 33, 34 |

Other column counts are rejected. The original format supports local source-data
reproduction. Its additional sensor block and identifiers are not public data.

## Reference files

Pulse CSVs require `time_s` and `GR_Pulse`. Pulse edges mark OUT/BACK passes.
GaitRite workbooks use the `Trial_NN` filename pattern and original export labels,
including `Left/Right Foot`, `Heel On`, `Toe Off ` (trailing space),
`FootFall Object #`, `Stride Length` and `Stride Time`.
Spatial measurements are converted from centimetres to metres.

## Primary outputs

| File | Contents |
| --- | --- |
| stride_level_pairs.csv | Accepted stride pairs, study code, condition, foot, trial, pass and reference/estimated measures |
| alignment_qc.csv | Per-pass alignment decisions and match quality |
| IMU_Synced_results.json | Per-record processing summaries |
| IMU_Synced_accuracy_summary.csv | Event accuracy by condition |
| IMU_Synced_per_subject_detection.csv | Participant/condition detection summaries |
| IMU_Synced_spatiotemporal_summary.csv | Timing and cadence summaries |
| IMU_Synced_stride_length_summary.csv | Estimated stride-length summaries |

Pair-table distances are metres and durations are seconds. Participant-level
outputs are ignored by Git and are not included in release assets.
