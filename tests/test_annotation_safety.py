"""Fail closed on truncated recipes and pipeline errors, without rescanning BCF."""

import logging
import subprocess
from unittest.mock import Mock

import pytest

from vcfcache.database.workflow_manager import BcftoolsCommand, WorkflowManager


@pytest.mark.parametrize(
    "command", ["false | true", "false\ntrue", "echo $UNDEFINED_SAFETY_TEST_VARIABLE"]
)
def test_actual_shell_enforces_saved_script(tmp_path, command):
    work = tmp_path / "path with spaces"
    work.mkdir()
    with pytest.raises(subprocess.CalledProcessError):
        BcftoolsCommand(command, logging.getLogger(__name__), work).run()
    assert (work / "stderr.txt").exists()
    assert "set -euo pipefail" in (work / "command.sh").read_text()


def test_unchecked_shell_still_reports_failure(tmp_path):
    result = BcftoolsCommand("false | true", logging.getLogger(__name__), tmp_path).run(
        check=False
    )
    assert result.returncode != 0


def manager(counts, allow=False):
    obj = WorkflowManager.__new__(WorkflowManager)
    obj.params_file_content = {
        "bcftools_cmd": "bcftools",
        "allow_annotation_drops": allow,
    }
    obj.logger = Mock()
    obj._count_variants = Mock(side_effect=counts)
    return obj


@pytest.mark.parametrize(
    "counts", [(100, 16), (100, 0), (100, 101), (100, None), (None, 100)]
)
def test_annotation_count_fails_closed(tmp_path, counts):
    with pytest.raises(RuntimeError):
        manager(counts)._validate_annotation_count(
            tmp_path / "input", tmp_path / "output"
        )


def test_intentional_filtering_requires_explicit_boolean(tmp_path):
    obj = manager([100, 16], allow=True)
    assert obj._validate_annotation_count(tmp_path / "in", tmp_path / "out") == (
        100,
        16,
    )
    obj.logger.warning.assert_called_once()
    with pytest.raises(ValueError, match="boolean"):
        manager([100, 16], allow="false")._validate_annotation_count(
            tmp_path / "in", tmp_path / "out"
        )
    with pytest.raises(RuntimeError):
        manager([100, None], allow=True)._validate_annotation_count(
            tmp_path / "in", tmp_path / "out"
        )


def test_complete_output_and_empty_input(tmp_path):
    for count in (0, 100):
        assert manager([count, count])._validate_annotation_count(
            tmp_path / "in", tmp_path / "out"
        ) == (count, count)
