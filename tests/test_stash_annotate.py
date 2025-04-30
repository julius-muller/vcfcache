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



def test_stash_annotate_with_add(test_output_dir):
    """Test the stash-annotate workflow with stash-add."""
    print("\n=== Testing stash-annotate workflow with stash-add ===")

    # Step 1: Run stash-init
    print("Running stash-init...")
    init_result = run_stash_init(TEST_VCF, test_output_dir, force=True)
    assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"

    # Step 2: Run stash-add
    print("Running stash-add...")
    add_result = run_stash_add(test_output_dir, TEST_VCF2)
    assert add_result.returncode == 0, f"stash-add failed: {add_result.stderr}"

    # Print information about the workflow directory and files
    workflow_dir = os.path.join(test_output_dir, "workflow")
    print(f"Workflow directory exists: {os.path.exists(workflow_dir)}")
    if os.path.exists(workflow_dir):
        print(f"Workflow directory contents: {os.listdir(workflow_dir)}")

    # Step 3: Run stash-annotate
    print("Running stash-annotate...")
    annotate_name = "test_annotation"
    annotate_result = run_stash_annotate(test_output_dir, annotate_name, force=True)
    if annotate_result.returncode != 0:
        print(f"Command output: {annotate_result.stdout}")
        print(f"Command error: {annotate_result.stderr}")
        print(f"Working directory contents: {os.listdir(test_output_dir)}")
        print(f"Workflow directory contents: {os.listdir(workflow_dir)}")
    assert annotate_result.returncode == 0, f"stash-annotate failed: {annotate_result.stderr}"

    # Get bcftools path for verification
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Step 4: Verify the annotation directory was created
    stash_dir = os.path.join(test_output_dir, "stash")
    annotation_dir = os.path.join(stash_dir, annotate_name)
    assert os.path.exists(annotation_dir), f"Annotation directory not found: {annotation_dir}"
    print(f"Annotation directory created: {annotation_dir}")

    # Step 5: Verify the annotated BCF file was created
    annotated_file = os.path.join(annotation_dir, "vcfstash_annotated.bcf")
    assert os.path.exists(annotated_file), f"Annotated BCF file not found: {annotated_file}"
    print(f"Annotated BCF file created: {annotated_file}")

    # Step 6: Verify the annotation.config file was copied
    config_file = os.path.join(annotation_dir, "annotation.config")
    assert os.path.exists(config_file), f"Annotation config file not found: {config_file}"
    print(f"Annotation config file created: {config_file}")

    # Step 7: Verify the blueprint_snapshot.info file was created
    snapshot_file = os.path.join(annotation_dir, "blueprint_snapshot.info")
    assert os.path.exists(snapshot_file), f"Blueprint snapshot file not found: {snapshot_file}"
    print(f"Blueprint snapshot file created: {snapshot_file}")

    # Step 8: Verify the annotated BCF file has the MOCK_ANNO tag
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", annotated_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Failed to view annotated file: {view_result.stderr}"
    assert "MOCK_ANNO" in view_result.stdout, "MOCK_ANNO tag not found in the header"
    print("MOCK_ANNO tag found in header")

    print("Successfully tested stash-annotate workflow with stash-add")
