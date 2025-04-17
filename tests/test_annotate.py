"""Test annotation functionality of VCFstash."""

import os
import sys
import pytest
from pathlib import Path
import subprocess
import shutil
import tempfile

# Constants
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TEST_VCF = os.path.join(TEST_DATA_DIR, "nodata", "crayz_db.bcf")
TEST_CONFIG = os.path.join(os.path.dirname(__file__), "config", "nextflow_test.config")
TEST_ANNO_CONFIG = os.path.join(os.path.dirname(__file__), "config", "annotation.config")
VCFSTASH_CMD = os.path.join(os.path.dirname(os.path.dirname(__file__)), "vcfstash.py")
EXPECTED_OUTPUT_DIR = os.path.join(TEST_DATA_DIR, "expected_output")

def test_cached_vs_uncached_annotation():
    """Test that cached and uncached annotation results match except for headers."""
    cached_bcf = Path(EXPECTED_OUTPUT_DIR) / "annotate_result" / "cached" / "sample4_vst.bcf"
    uncached_bcf = Path(EXPECTED_OUTPUT_DIR) / "annotate_result" / "uncached" / "sample4_vst.bcf"

    assert cached_bcf.exists(), f"Cached reference BCF not found at {cached_bcf}"
    assert uncached_bcf.exists(), f"Uncached reference BCF not found at {uncached_bcf}"

    # Read BCF contents excluding headers
    def get_variants(bcf_path):
        try:
            # Get bcftools path
            bcftools_path = os.path.join(os.environ.get('VCFSTASH_ROOT', ''), 'tools', 'bcftools')
            if not os.path.exists(bcftools_path):
                # Fall back to system bcftools if the project-specific one doesn't exist
                bcftools_path = 'bcftools'

            bcf_text = subprocess.run(
                [bcftools_path, "view", str(bcf_path)],
                capture_output=True,
                text=True,
                check=True
            ).stdout
            return [line for line in bcf_text.splitlines()
                    if line and not line.startswith('#')]
        except subprocess.CalledProcessError as e:
            pytest.fail(f"Failed to read BCF file {bcf_path}: {e.stderr}")
            return []

    cached_variants = get_variants(cached_bcf)
    uncached_variants = get_variants(uncached_bcf)

    assert len(cached_variants) > 0, f"No variants found in cached BCF {cached_bcf}"
    assert len(uncached_variants) > 0, f"No variants found in uncached BCF {uncached_bcf}"

    # Find and display the first difference
    for i, (cached, uncached) in enumerate(zip(cached_variants, uncached_variants)):
        if cached != uncached:
            pytest.fail(
                f"First difference at line {i + 1}:\n"
                f"Cached:   {cached}\n"
                f"Uncached: {uncached}"
            )

    # If lengths differ, show which file has extra lines
    if len(cached_variants) != len(uncached_variants):
        pytest.fail(
            f"Number of variants differs: "
            f"cached={len(cached_variants)}, uncached={len(uncached_variants)}"
        )

@pytest.mark.parametrize("use_cache", [True, False])
def test_annotate_command(use_cache):
    """Test that the annotate command works correctly with and without cache."""
    # Skip if reference files don't exist
    test_input = Path(TEST_DATA_DIR) / "sample4.bcf"
    test_cache = Path(EXPECTED_OUTPUT_DIR) / "stash_add_annotate_result/stash/test_annotation"

    if not test_input.exists() or not test_cache.exists():
        pytest.skip(f"Required test files not found: {test_input} or {test_cache}")

    # Create temporary output directory
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "output"
        output_dir.mkdir()

        # Build command
        cmd = [
            sys.executable, VCFSTASH_CMD,
            "annotate",
            "-i", str(test_input),
            "-a", str(test_cache),
            "-o", str(output_dir),
            "-f"
        ]

        # Add uncached flag if testing without cache
        if not use_cache:
            cmd.append("--uncached")

        # Run command
        try:
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                check=True,
                timeout=300
            )

            # Check that output file exists
            output_file = output_dir / "sample4_vst.bcf"
            assert output_file.exists(), f"Output file not created: {output_file}"

            # Check that file has content
            file_size = output_file.stat().st_size
            assert file_size > 0, f"Output file is empty: {output_file}"

        except subprocess.CalledProcessError as e:
            pytest.fail(f"Command failed with exit code {e.returncode}:\n{e.stderr}")
        except subprocess.TimeoutExpired:
            pytest.fail("Command timed out after 5 minutes")
