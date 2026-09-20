# Algorithm and code map

Start with `scripts/run_pipeline.py`. It coordinates one participant-condition
record at a time. Functions in `src/` perform the signal processing. The primary
numerical engine is preserved from submitted commit `7530ddc`.

| Stage | Implementation | Purpose |
| --- | --- | --- |
| Load | `run_pipeline.load_imu_csv`, `config.imu_column_layout` | Read 29- or 43-column exports, infer sampling rate and filter sagittal gyroscope signals. |
| Detect | `gait_events/adaptive.py` | Detect heel strikes and toe-offs using amplitude and cadence information. |
| Integrate | `stride_length.estimate_stride_length_zupt` | Confirm stance, estimate orientation, remove gravity and integrate with zero-velocity and residual-gravity corrections. |
| Synchronize | `time_sync.load_gr_pulse` | Map pulse timing to the IMU clock. |
| Align | `time_sync.find_optimal_alignment` | Require at least 60% heel-strike matches within 150 ms for each pass. |
| Pair | `run_pipeline.py` | Match same-foot, on-mat strides by reciprocal nearest midpoints and common plausibility bounds. |
| Summarize | `run_stats.py`, `validate_release.py` | Compute agreement and check frozen reference targets. |

## Interpretation

- Integration uses mid-stance-to-mid-stance intervals. Reported heel-strike
  anchors support matching, and are not the displacement integration boundaries.
- Only the one-dimensional left gyroscope signal is inverted for event detection.
  Three-dimensional frames are preserved for orientation estimation.
- Failed alignment contributes no paired strides. The documented threshold-only
  stance fallback is distinct from alignment failure.
- Fixed named walk exclusions in `src/config.py` address documented data-quality
  issues. No participant-condition exclusion depends on error or pair count.
- `run_contralateral_vs_unilateral.py` is a supplementary stance comparison,
  with four-decimal pair outputs. Use the primary pipeline's six-decimal table
  to reproduce headline results.

The validation covers accepted supervised straight passes. It does not establish
accuracy during turns, everyday walking, different placements or other cohorts.
Changing thresholds or integration settings creates a new analysis.
