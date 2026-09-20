"""Synthetic interface checks for the public Paper 1 code release."""

import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from comparison import match_events  # noqa: E402
from config import imu_column_layout, subject_code  # noqa: E402


def load_validator():
    spec = importlib.util.spec_from_file_location(
        "validate_release", ROOT / "scripts" / "validate_release.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("column_count", [29, 43])
def test_supported_imu_schemas_include_bilateral_sensor_channels(column_count):
    layout = imu_column_layout(column_count)
    assert len(layout["left_accel"]) == len(layout["right_accel"]) == 3
    assert len(layout["left_gyro"]) == len(layout["right_gyro"]) == 3
    assert layout["left_gyro_z"] in layout["left_gyro"]
    assert layout["right_gyro_z"] in layout["right_gyro"]


def test_unsupported_imu_schema_names_supported_counts():
    with pytest.raises(ValueError, match="Expected 43 or 29"):
        imu_column_layout(42)


def test_subject_code_uses_only_final_raw_folder_token():
    assert subject_code("EXAMPLE_01_01_00_123") == "ID_123"
    assert subject_code("ID_107") == "ID_107"


def test_event_matching_is_unique_and_respects_tolerance():
    matched, missed, false_pos = match_events(
        np.array([0.98, 1.03, 2.30]), np.array([1.00, 2.00]), tolerance_s=0.05
    )
    assert matched == [(1.0, 0.98)]
    assert missed == [2.0]
    assert false_pos == [1.03, 2.3]


@pytest.mark.parametrize(
    "script_name",
    ["run_pipeline.py", "run_contralateral_vs_unilateral.py"],
)
def test_cli_rejects_unknown_condition_before_data_loading(script_name):
    result = subprocess.run(
        [
            sys.executable,
            ROOT / "scripts" / script_name,
            "--base-dir",
            "/missing",
            "--conditions",
            "BAD",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 2
    assert "unknown condition code(s): BAD" in result.stderr


@pytest.mark.parametrize(
    "script_name",
    ["run_pipeline.py", "run_contralateral_vs_unilateral.py"],
)
def test_cli_rejects_missing_data_root(script_name, tmp_path):
    missing_root = tmp_path / "missing"
    result = subprocess.run(
        [sys.executable, ROOT / "scripts" / script_name, "--base-dir", missing_root],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 2
    assert f"data root does not exist: {missing_root}" in result.stderr


@pytest.mark.parametrize(
    "script_name",
    ["run_pipeline.py", "run_contralateral_vs_unilateral.py"],
)
def test_cli_rejects_data_root_with_missing_subdirectories(script_name, tmp_path):
    (tmp_path / "IMU Data").mkdir()
    result = subprocess.run(
        [sys.executable, ROOT / "scripts" / script_name, "--base-dir", tmp_path],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 2
    assert "data root is missing required subdirectories" in result.stderr


@pytest.mark.parametrize(
    "script_name",
    ["run_pipeline.py", "run_contralateral_vs_unilateral.py"],
)
def test_cli_rejects_root_with_no_usable_subjects(script_name, tmp_path):
    for directory in ["IMU Data", "GaitRite Data", "Synchronization Data"]:
        (tmp_path / directory).mkdir()
    result = subprocess.run(
        [sys.executable, ROOT / "scripts" / script_name, "--base-dir", tmp_path],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 2
    assert "no usable subject directories found under IMU Data" in result.stderr


@pytest.mark.parametrize(
    "script_name",
    ["run_pipeline.py", "run_contralateral_vs_unilateral.py"],
)
def test_cli_rejects_requested_subject_that_is_not_present(script_name, tmp_path):
    for directory in ["IMU Data", "GaitRite Data", "Synchronization Data"]:
        (tmp_path / directory).mkdir()
    result = subprocess.run(
        [
            sys.executable,
            ROOT / "scripts" / script_name,
            "--base-dir",
            tmp_path,
            "--subjects",
            "missing_subject",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 2
    assert "subject(s) not found under IMU Data: missing_subject" in result.stderr


def test_validator_rejects_non_boolean_alignment_values():
    validator = load_validator()
    with pytest.raises(validator.ReleaseValidationError, match="only true or false"):
        validator.parse_valid_column(pd.Series(["true", "yes"]))


def test_validator_rejects_wrong_counts_even_with_python_optimization(tmp_path):
    pd.DataFrame(columns=["condition", "subject", "zupt_sl_m", "gr_sl_m"]).to_csv(
        tmp_path / "stride_level_pairs.csv", index=False
    )
    pd.DataFrame(columns=["condition", "valid"]).to_csv(
        tmp_path / "alignment_qc.csv", index=False
    )
    result = subprocess.run(
        [
            sys.executable,
            "-O",
            ROOT / "scripts" / "validate_release.py",
            "--data-dir",
            tmp_path,
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode != 0
    assert "stride pair count: observed 0" in result.stderr
