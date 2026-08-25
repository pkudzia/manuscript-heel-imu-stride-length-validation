"""Constants shared by the synced-data pipelines (run_pipeline.py and
run_contralateral_vs_unilateral.py): sensor layout, filter settings,
GaitRite plausibility bounds, and the data-quality exclusion lists."""


def subject_code(folder_name):
    """De-identified study code for a raw-data folder name.

    Raw data folders end in the numeric study code. All pipeline outputs
    and the exclusion lists below use this ID_NNN form, which matches the
    participant labels in the manuscript and supplement.
    """
    return "ID_" + folder_name.rsplit("_", 1)[-1]


# Do not exclude participant-condition records because their stride-length
# error is large or because they contribute few pairs. Those are study
# outcomes. Conditions without the required source files are skipped during
# loading and are reported separately.
EXCLUDE_STRIDE_PAIRS = set()

# IMU sensor settings (APDM Opal, 256 Hz, heel-mounted bilaterally)
SAMPLING_FREQ = 256
BANDPASS_LOW = 0.5
BANDPASS_HIGH = 10.0
FILTER_ORDER = 4

# IMU column indices
RIGHT_GYRO_Z_COL = 6
LEFT_GYRO_Z_COL = 34
RIGHT_ACCEL_COLS = [1, 2, 3]
RIGHT_GYRO_COLS = [4, 5, 6]
LEFT_ACCEL_COLS = [29, 30, 31]
LEFT_GYRO_COLS = [32, 33, 34]

# Polarity convention: invert left foot Z gyro so both feet share a
# common sign convention for heel-strike detection. Applied ONLY to the
# 1D filtered Z signal used by detectors. The 3-axis arrays consumed by
# Mahony AHRS / ZUPT preserve the raw sensor frame (a single-axis flip
# would create a left-handed frame and corrupt orientation tracking).
INVERT_LEFT_GYRO = True
INVERT_RIGHT_GYRO = False

# GR stride length plausibility bounds (m). Values outside this range
# are GaitRite measurement errors (e.g. missed footfalls creating
# double-length strides) and are excluded from pairing.
GR_SL_MIN_M = 0.30
GR_SL_MAX_M = 2.00

# GR stride TIME plausibility bounds (s). Adult walking stride time is
# 1.0-1.2 s normally, up to ~1.6-1.7 s for slow / shuffling MCI.
# Values outside this range are GaitRite measurement errors — typically
# a missed intermediate heel strike that causes GR to report two real
# strides as a single "double stride" of ~2x the true duration.
# Filter rejects these so we don't pair a known-broken GR reference
# value to a correctly-detected IMU stride.
GR_ST_MIN_S = 0.40
GR_ST_MAX_S = 1.80

# Walk-level exclusions: specific (subject_code, condition, trial_num,
# walk) tuples to drop from pairing. Use sparingly and document each
# entry. Each is a walk where a real IMU-vs-GR comparison is not
# possible due to data-quality issues that are not the algorithm's
# responsibility.
EXCLUDE_WALKS = {
    # ID_107 UL trial 2 OUT: IMU heel-strike detector picked up a
    # pre-mat heel strike (~50.7 s, before the participant entered the
    # mat). ZUPT then integrated across the mat-entry transition,
    # producing a ~0.8 m "stride" that does not correspond to any GR
    # stride.
    ("ID_107", "UL", 2, "OUT"),
    # ID_149 L2 trial 1 BACK: same mat-entry pattern. IMU detected a
    # pre-mat heel strike (~74.9 s); GR captured only on-mat strides
    # starting with a short gait-initiation stride (0.98 m). The
    # mismatch produced a +47 cm error.
    ("ID_149", "L2", 1, "BACK"),
    # ID_195 L7 trial 1 BACK: GaitRite internal data inconsistency. The
    # exported HS times show only 2 events 2.0 s apart for the right
    # foot, but the spatiotemporal export reports stride_time = 1.475 s
    # for the same stride. GR's HS detection and stride-time
    # computation disagree, so any pairing here is unreliable.
    ("ID_195", "L7", 1, "BACK"),
}

# Reported conditions: condition code -> (folder name, GaitRite filename label)
CONDITIONS = {
    "UL": ("Unloaded", "Unloaded"),
    "L2": ("Loaded 2", "Loaded2"),
    "L3": ("Loaded 3", "Loaded3"),
    "L7": ("Loaded 7", "Loaded7"),
}
