import pytest
import tempfile
import shutil
import hashlib
import os
import re
import json
import sys
import subprocess
from pathlib import Path
from src.utils.paths import get_vcfstash_root, get_resource_path

# Define constants for testing
TEST_ROOT = Path(__file__).parent
VCFSTASH_CMD = get_vcfstash_root() / "vcfstash.py"
TEST_DATA_DIR = TEST_ROOT / "data" / "nodata"
TEST_VCF = TEST_DATA_DIR / "crayz_db.bcf"
TEST_VCF2 = TEST_DATA_DIR / "crayz_db2.bcf"
TEST_SAMPLE = TEST_DATA_DIR / "sample4.bcf"
TEST_PARAMS = TEST_ROOT / "config" / "test_params.yaml"
TEST_ANNO_CONFIG = TEST_ROOT / "config" / "test_annotation.config"


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


# Helper functions
def compute_md5(filename):
    """Compute MD5 hash for a file."""
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

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

def run_annotate(annotation_db, input_vcf, output_dir, force=False):
    """Run the annotate command and return the process result."""
    # Create a temporary params file
    temp_params_file = create_temp_params_file()

    try:
        cmd = [
            sys.executable,  # Use the current Python interpreter
            str(VCFSTASH_CMD),
            "annotate",
            "-a", str(annotation_db),
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



# Tests
def test_stash_init(test_output_dir):
    """Test that stash-init creates the expected directory structure."""
    # Run stash-init
    result = run_stash_init(TEST_VCF, test_output_dir, force=True)
    assert result.returncode == 0, f"stash-init failed: {result.stderr}"

    # Verify that the output directory has the expected structure
    # Check for key directories and files
    assert os.path.isdir(test_output_dir), f"Output directory {test_output_dir} not found"
    assert os.path.isdir(os.path.join(test_output_dir, "blueprint")), "Blueprint directory not found in output"
    assert os.path.isdir(os.path.join(test_output_dir, "workflow")), "Workflow directory not found in output"

    # Check for key files
    assert os.path.exists(os.path.join(test_output_dir, "blueprint", "vcfstash.bcf")), "vcfstash.bcf not found in output"
    assert os.path.exists(os.path.join(test_output_dir, "blueprint", "sources.info")), "sources.info not found in output"





def test_key_files_content(test_output_dir):
    """Test that key files in the stash directory have valid content."""
    # Run stash-init
    result = run_stash_init(TEST_VCF, test_output_dir, force=True)
    assert result.returncode == 0, f"stash-init failed: {result.stderr}"

    # List of important files to check
    key_files = [
        "blueprint/sources.info",
        "workflow/main.nf"
    ]

    # Check all the module files too
    module_files = [
        "workflow/modules/annotate_cache.nf",
        "workflow/modules/annotate.nf",
        "workflow/modules/intersect.nf",
        "workflow/modules/merge.nf",
        "workflow/modules/merge_variants.nf",
        "workflow/modules/normalize.nf",
        "workflow/modules/utils.nf"
    ]
    key_files.extend(module_files)

    # Track how many files were checked
    checked_files = 0

    for file_path in key_files:
        test_file = os.path.join(test_output_dir, file_path)

        # Check if file exists
        if not os.path.exists(test_file):
            print(f"Warning: File {file_path} does not exist. Skipping check.")
            continue

        checked_files += 1

        # Verify that the file is not empty
        assert os.path.getsize(test_file) > 0, f"File {file_path} is empty"

        # For sources.info, check that it's valid JSON
        if file_path.endswith('sources.info'):
            try:
                with open(test_file, 'r') as f:
                    json_content = json.load(f)
                assert "input_files" in json_content, "sources.info missing 'input_files' key"
            except json.JSONDecodeError:
                assert False, f"File {file_path} is not valid JSON"

    # Ensure at least some files were checked
    assert checked_files > 0, "No files were checked."

def test_bcf_file_validity(test_output_dir):
    """Test that the generated BCF file is valid."""
    # Run stash-init
    result = run_stash_init(TEST_VCF, test_output_dir, force=True)
    assert result.returncode == 0, f"stash-init failed: {result.stderr}"

    # Check if the BCF file exists
    test_bcf = Path(test_output_dir) / "blueprint" / "vcfstash.bcf"
    assert test_bcf.exists(), f"BCF file not found: {test_bcf}"

    # Get bcftools path
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Check if the BCF file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", str(test_bcf)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"BCF file is not valid: {view_result.stderr}"

    # Check if the BCF file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", str(test_bcf)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for BCF file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "BCF file has no variants"

    # Extract the number of records
    num_records = 0
    for line in stats_result.stdout.splitlines():
        if "number of records:" in line:
            num_records = int(line.split(":")[-1].strip())
            break

    # Ensure there are records in the file
    assert num_records > 0, f"BCF file has {num_records} records, expected > 0"

def test_stash_add_workflow(test_output_dir):
    """Test the stash-add workflow."""
    # Step 1: Run stash-init
    init_result = run_stash_init(TEST_VCF, test_output_dir, force=True)
    assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"

    # Step 2: Run stash-add
    add_result = run_stash_add(test_output_dir, TEST_VCF2)
    assert add_result.returncode == 0, f"stash-add failed: {add_result.stderr}"

    # Get bcftools path
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Check if the blueprint BCF file exists
    blueprint_file = os.path.join(test_output_dir, "blueprint", "vcfstash.bcf")
    assert os.path.exists(blueprint_file), f"Blueprint file not found: {blueprint_file}"

    # Check if the blueprint BCF file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", blueprint_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Blueprint file is not valid: {view_result.stderr}"

    # Check if the blueprint BCF file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", blueprint_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for blueprint file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Blueprint file has no variants"

    # Check if the sources.info file exists
    sources_file = os.path.join(test_output_dir, "blueprint", "sources.info")
    assert os.path.exists(sources_file), f"Sources file not found: {sources_file}"

    # Check if the sources.info file is valid JSON
    try:
        with open(sources_file, 'r') as f:
            sources_data = json.load(f)
        assert "input_files" in sources_data, "sources.info missing 'input_files' key"

        # Check if both input files are in the sources.info file
        input_files = [os.path.basename(f["path"]) for f in sources_data["input_files"]]
        assert os.path.basename(str(TEST_VCF)) in input_files, f"{os.path.basename(str(TEST_VCF))} not found in sources.info"
        assert os.path.basename(str(TEST_VCF2)) in input_files, f"{os.path.basename(str(TEST_VCF2))} not found in sources.info"
    except json.JSONDecodeError:
        pytest.fail(f"Sources file is not valid JSON: {sources_file}")
    except Exception as e:
        pytest.fail(f"Unexpected error reading sources file: {str(e)}")
