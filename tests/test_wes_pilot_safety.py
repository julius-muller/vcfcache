"""The pilot actually executes the recorded strict shell script."""

import logging
import subprocess

import pytest

from vcfcache.database.workflow_manager import BcftoolsCommand


@pytest.mark.parametrize("script", ["false | true", "false\ntrue"])
def test_pipeline_errors_cannot_be_masked(tmp_path, script):
    directory = tmp_path / "directory with spaces"
    directory.mkdir()
    with pytest.raises(subprocess.CalledProcessError):
        BcftoolsCommand(script, logging.getLogger(__name__), directory).run()
    assert (directory / "stderr.txt").exists()
