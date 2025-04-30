"""Test stash-annotate functionality of VCFstash."""

import os
import sys
import json
import pytest
from pathlib import Path
import subprocess
import shutil
import tempfile
from src.utils.paths import get_vcfstash_root, get_resource_path

# Constants
TEST_ROOT = Path(__file__).parent
TEST_DATA_DIR = TEST_ROOT / "data" / "nodata"
TEST_VCF = TEST_DATA_DIR / "crayz_db.bcf"
TEST_VCF2 = TEST_DATA_DIR / "crayz_db2.bcf"
TEST_SAMPLE = TEST_DATA_DIR / "sample4.bcf"
TEST_PARAMS = TEST_ROOT / "config" / "test_params.yaml"
TEST_ANNO_CONFIG = TEST_ROOT / "config" / "test_annotation.config"
VCFSTASH_CMD = get_vcfstash_root() / "vcfstash.py"


# Fixture for temporary test output directory
@pytest.fixture
def test_output_dir():
    """Provides a path for a directory that doesn't exist yet."""
    # Create a path for a temporary directory, but don't create it
    temp_dir = tempfile.mkdtemp(prefix="vcfstash_test_")

    # Remove the directory immediately - vcfstash will create it
    os.rmdir(temp_dir)

    # Return the path to the test function
    yield temp_dir

    # Clean up after the test is done
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)


def create_temp_params_file():
    """Create a temporary params file with the correct paths and return its path."""
    # Get the VCFSTASH_ROOT directory
    vcfstash_root = str(get_vcfstash_root())

    # Create a temporary params file with the correct paths
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
    temp_params_file = temp_file.name

    # Read the original params file
    with open(TEST_PARAMS, 'r') as f:
        params_content = f.read()

    # Replace ${VCFSTASH_ROOT} with the actual value
    params_content = params_content.replace('${VCFSTASH_ROOT}', vcfstash_root)

    # Write the modified content to the temporary file
    temp_file.write(params_content)
    temp_file.close()

    return temp_params_file


def run_stash_init(input_vcf, output_dir, force=False):
    """Run the stash-init command and return the process result."""
    # Make sure the directory doesn't exist (clean start)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    # Create a temporary params file
    temp_params_file = create_temp_params_file()

    try:
        cmd = [
            sys.executable,  # Use the current Python interpreter
            str(VCFSTASH_CMD),
            "stash-init",
            "--vcf", str(input_vcf),
            "--output", str(output_dir),
            "-y", temp_params_file
        ]

        if force:
            cmd.append("-f")

        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result
    finally:
        # Clean up the temporary file
        if os.path.exists(temp_params_file):
            os.unlink(temp_params_file)


def run_stash_add(db_dir, input_vcf):
    """Run the stash-add command and return the process result."""
    cmd = [
        sys.executable,  # Use the current Python interpreter
        str(VCFSTASH_CMD),
        "stash-add",
        "--db", str(db_dir),
        "-i", str(input_vcf)
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return result


def run_stash_annotate(db_dir, name, force=False):
    """Run the stash-annotate command and return the process result."""
    # Create a temporary params file
    temp_params_file = create_temp_params_file()

    try:
        cmd = [
            sys.executable,  # Use the current Python interpreter
            str(VCFSTASH_CMD),
            "stash-annotate",
            "--name", name,
            "-a", str(TEST_ANNO_CONFIG),
            "--db", str(db_dir),
            "-y", temp_params_file
        ]

        if force:
            cmd.append("-f")

        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result
    finally:
        # Clean up the temporary file
        if os.path.exists(temp_params_file):
            os.unlink(temp_params_file)


def test_stash_annotate_workflow(test_output_dir):
    """Test the stash-annotate workflow."""
    # Skip this test for now as it requires a working Nextflow workflow
    pytest.skip("Skipping stash-annotate test as it requires a working Nextflow workflow. "
                "The annotate functionality is tested directly in test_annotate.py.")


def test_stash_annotate_with_add(test_output_dir):
    """Test the stash-annotate workflow with stash-add."""
    # Skip this test for now as it requires a working Nextflow workflow
    pytest.skip("Skipping stash-annotate test as it requires a working Nextflow workflow. "
                "The annotate functionality with stash-add is tested directly in test_annotate.py.")
