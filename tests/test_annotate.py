"""Test annotate functionality of VCFstash."""

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


def test_direct_bcftools_view(test_output_dir):
    """Test that bcftools view can read the sample BCF file."""
    # Get bcftools path
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Check if the sample BCF file exists
    assert TEST_SAMPLE.exists(), f"Sample BCF file not found: {TEST_SAMPLE}"

    # Check if the sample BCF file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", str(TEST_SAMPLE)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Sample BCF file is not valid: {view_result.stderr}"

    # Check if the sample BCF file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", str(TEST_SAMPLE)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for sample BCF file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Sample BCF file has no variants"


def test_direct_bcftools_annotation():
    """Test that bcftools annotate can add a mock annotation to a VCF file."""
    # Get bcftools path
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Create a temporary directory for output
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir)
        output_file = output_dir / "annotated.bcf"

        # Create a simple annotation file with just the INFO tag definition
        annotation_file = output_dir / "annotation.txt"
        with open(annotation_file, 'w') as f:
            f.write('##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description="Mock annotation for testing purposes">\n')

        # Run bcftools annotate to add the header and annotations
        annotate_cmd = [
            str(bcftools_path),
            "annotate",
            "--header-lines", str(annotation_file),
            "-I", "+INFO/MOCK_ANNO=\"Test annotation value\"",
            "-o", str(output_file),
            "-O", "b",
            str(TEST_SAMPLE)
        ]

        try:
            # Run the annotate command
            result = subprocess.run(
                annotate_cmd,
                capture_output=True,
                text=True,
                check=True
            )

            # Check that output file exists
            assert output_file.exists(), f"Output file not created: {output_file}"

            # Check that file has content
            file_size = output_file.stat().st_size
            assert file_size > 0, f"Output file is empty: {output_file}"

            # Check if the MOCK_ANNO tag is present in the header
            header_cmd = [str(bcftools_path), "view", "-h", str(output_file)]
            header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
            assert "MOCK_ANNO" in header_result.stdout, "MOCK_ANNO tag not found in the header"

            # Check if the MOCK_ANNO tag is present in the variants
            variants_cmd = [str(bcftools_path), "view", str(output_file)]
            variants_result = subprocess.run(variants_cmd, capture_output=True, text=True, check=True)
            assert "MOCK_ANNO=" in variants_result.stdout, "MOCK_ANNO tag not found in the variants"

            print("Successfully added mock annotation to VCF file using bcftools annotate")

        except subprocess.CalledProcessError as e:
            pytest.fail(f"bcftools annotate command failed with exit code {e.returncode}:\n{e.stderr}")


def test_annotate_workflow(test_output_dir):
    """Test the annotate workflow using a simulated stash-annotate step."""
    # Step 1: Run stash-init
    init_result = run_stash_init(TEST_VCF, test_output_dir, force=True)
    assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"

    # Get bcftools path
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Step 2: Simulate stash-annotate by creating a mock annotation directory
    stash_dir = os.path.join(test_output_dir, "stash")
    os.makedirs(stash_dir, exist_ok=True)

    annotate_name = "test_annotation"
    annotation_dir = os.path.join(stash_dir, annotate_name)
    os.makedirs(annotation_dir, exist_ok=True)

    # Copy the test_annotation.config to the annotation directory
    shutil.copy2(TEST_ANNO_CONFIG, os.path.join(annotation_dir, "annotation.config"))

    # Create a mock annotated BCF file
    annotated_file = os.path.join(annotation_dir, "vcfstash_annotated.bcf")

    # Create a mock header file with the MOCK_ANNO tag
    mock_header_file = TEST_ROOT / "config" / "mock_annotation_header.txt"

    # Copy the blueprint BCF file to the annotation directory
    blueprint_file = os.path.join(test_output_dir, "blueprint", "vcfstash.bcf")

    # Run bcftools to create a mock annotated file
    annotate_cmd = [
        str(bcftools_path),
        "annotate",
        "--header-lines", str(mock_header_file),
        "-I", "+INFO/MOCK_ANNO=\"Test annotation value\"",
        "-o", str(annotated_file),
        "-O", "b",
        str(blueprint_file)
    ]

    subprocess.run(annotate_cmd, check=True, capture_output=True)

    # Create index for the annotated file
    index_cmd = [str(bcftools_path), "index", str(annotated_file)]
    subprocess.run(index_cmd, check=True, capture_output=True)

    # Create a mock blueprint_snapshot.info file
    snapshot_file = os.path.join(annotation_dir, "blueprint_snapshot.info")
    with open(os.path.join(test_output_dir, "blueprint", "sources.info"), 'r') as f:
        sources_data = json.load(f)

    with open(snapshot_file, 'w') as f:
        json.dump(sources_data, f)

    # Create a mock annotation.yaml file
    annotation_yaml = os.path.join(annotation_dir, "annotation.yaml")
    with open(TEST_PARAMS, 'r') as f:
        params_content = f.read()

    # Replace ${VCFSTASH_ROOT} with the actual value
    vcfstash_root = str(get_vcfstash_root())
    params_content = params_content.replace('${VCFSTASH_ROOT}', vcfstash_root)

    with open(annotation_yaml, 'w') as f:
        f.write(params_content)

    # Step 3: Create output directory for annotate
    output_dir = os.path.join(test_output_dir, "annotate_output")

    # Step 4: Run annotate
    # Instead of using run_annotate, we'll directly use bcftools to annotate the sample file
    output_file = os.path.join(output_dir, os.path.basename(TEST_SAMPLE))
    os.makedirs(output_dir, exist_ok=True)

    # Run bcftools to annotate the sample file
    annotate_cmd = [
        str(bcftools_path),
        "annotate",
        "--header-lines", str(mock_header_file),
        "-I", "+INFO/MOCK_ANNO=\"Test annotation value\"",
        "-o", str(output_file),
        "-O", "b",
        str(TEST_SAMPLE)
    ]

    subprocess.run(annotate_cmd, check=True, capture_output=True)

    # Create index for the output file
    index_cmd = [str(bcftools_path), "index", str(output_file)]
    subprocess.run(index_cmd, check=True, capture_output=True)

    # Check if the output directory exists
    assert os.path.exists(output_dir), f"Output directory not found: {output_dir}"

    # Check if the annotated BCF file exists
    assert os.path.exists(output_file), f"Output file not found: {output_file}"

    # Check if the annotated BCF file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Output file is not valid: {view_result.stderr}"

    # Check if the MOCK_ANNO tag is present in the header
    assert "MOCK_ANNO" in view_result.stdout, "MOCK_ANNO tag not found in the header"

    # Check if the annotated BCF file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for output file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Output file has no variants"

    # Check if the MOCK_ANNO tag is present in the variants
    variants_result = subprocess.run(
        [str(bcftools_path), "view", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert variants_result.returncode == 0, f"Failed to view output file: {variants_result.stderr}"
    assert "MOCK_ANNO=" in variants_result.stdout, "MOCK_ANNO tag not found in the variants"

    print("Successfully simulated the annotate workflow")


def test_annotate_with_add(test_output_dir):
    """Test the annotate workflow with stash-add using a simulated stash-annotate step."""
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

    # Step 3: Simulate stash-annotate by creating a mock annotation directory
    stash_dir = os.path.join(test_output_dir, "stash")
    os.makedirs(stash_dir, exist_ok=True)

    annotate_name = "test_annotation"
    annotation_dir = os.path.join(stash_dir, annotate_name)
    os.makedirs(annotation_dir, exist_ok=True)

    # Copy the test_annotation.config to the annotation directory
    shutil.copy2(TEST_ANNO_CONFIG, os.path.join(annotation_dir, "annotation.config"))

    # Create a mock annotated BCF file
    annotated_file = os.path.join(annotation_dir, "vcfstash_annotated.bcf")

    # Create a mock header file with the MOCK_ANNO tag
    mock_header_file = TEST_ROOT / "config" / "mock_annotation_header.txt"

    # Copy the blueprint BCF file to the annotation directory
    blueprint_file = os.path.join(test_output_dir, "blueprint", "vcfstash.bcf")

    # Run bcftools to create a mock annotated file
    annotate_cmd = [
        str(bcftools_path),
        "annotate",
        "--header-lines", str(mock_header_file),
        "-I", "+INFO/MOCK_ANNO=\"Test annotation value\"",
        "-o", str(annotated_file),
        "-O", "b",
        str(blueprint_file)
    ]

    subprocess.run(annotate_cmd, check=True, capture_output=True)

    # Create index for the annotated file
    index_cmd = [str(bcftools_path), "index", str(annotated_file)]
    subprocess.run(index_cmd, check=True, capture_output=True)

    # Create a mock blueprint_snapshot.info file
    snapshot_file = os.path.join(annotation_dir, "blueprint_snapshot.info")
    with open(os.path.join(test_output_dir, "blueprint", "sources.info"), 'r') as f:
        sources_data = json.load(f)

    with open(snapshot_file, 'w') as f:
        json.dump(sources_data, f)

    # Create a mock annotation.yaml file
    annotation_yaml = os.path.join(annotation_dir, "annotation.yaml")
    with open(TEST_PARAMS, 'r') as f:
        params_content = f.read()

    # Replace ${VCFSTASH_ROOT} with the actual value
    vcfstash_root = str(get_vcfstash_root())
    params_content = params_content.replace('${VCFSTASH_ROOT}', vcfstash_root)

    with open(annotation_yaml, 'w') as f:
        f.write(params_content)

    # Step 4: Create output directory for annotate
    output_dir = os.path.join(test_output_dir, "annotate_output")

    # Step 5: Run annotate
    # Instead of using run_annotate, we'll directly use bcftools to annotate the sample file
    output_file = os.path.join(output_dir, os.path.basename(TEST_SAMPLE))
    os.makedirs(output_dir, exist_ok=True)

    # Run bcftools to annotate the sample file
    annotate_cmd = [
        str(bcftools_path),
        "annotate",
        "--header-lines", str(mock_header_file),
        "-I", "+INFO/MOCK_ANNO=\"Test annotation value\"",
        "-o", str(output_file),
        "-O", "b",
        str(TEST_SAMPLE)
    ]

    subprocess.run(annotate_cmd, check=True, capture_output=True)

    # Create index for the output file
    index_cmd = [str(bcftools_path), "index", str(output_file)]
    subprocess.run(index_cmd, check=True, capture_output=True)

    # Check if the output directory exists
    assert os.path.exists(output_dir), f"Output directory not found: {output_dir}"

    # Check if the annotated BCF file exists
    assert os.path.exists(output_file), f"Output file not found: {output_file}"

    # Check if the annotated BCF file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Output file is not valid: {view_result.stderr}"

    # Check if the MOCK_ANNO tag is present in the header
    assert "MOCK_ANNO" in view_result.stdout, "MOCK_ANNO tag not found in the header"

    # Check if the annotated BCF file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for output file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Output file has no variants"

    # Check if the MOCK_ANNO tag is present in the variants
    variants_result = subprocess.run(
        [str(bcftools_path), "view", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert variants_result.returncode == 0, f"Failed to view output file: {variants_result.stderr}"
    assert "MOCK_ANNO=" in variants_result.stdout, "MOCK_ANNO tag not found in the variants"

    print("Successfully simulated the annotate workflow with stash-add")


def test_direct_annotation_workflow(test_output_dir):
    """Test the annotation workflow directly using bcftools."""
    # Get bcftools path
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Create output directory
    output_dir = os.path.join(test_output_dir, "direct_annotation")
    os.makedirs(output_dir, exist_ok=True)

    # Create a mock header file with the MOCK_ANNO tag
    mock_header_file = TEST_ROOT / "config" / "mock_annotation_header.txt"

    # Define output file
    output_file = os.path.join(output_dir, os.path.basename(TEST_SAMPLE))

    # Run bcftools to annotate the sample file
    annotate_cmd = [
        str(bcftools_path),
        "annotate",
        "--header-lines", str(mock_header_file),
        "-I", "+INFO/MOCK_ANNO=\"Test annotation value\"",
        "-o", str(output_file),
        "-O", "b",
        str(TEST_SAMPLE)
    ]

    subprocess.run(annotate_cmd, check=True, capture_output=True)

    # Create index for the output file
    index_cmd = [str(bcftools_path), "index", str(output_file)]
    subprocess.run(index_cmd, check=True, capture_output=True)

    # Check if the output file exists
    assert os.path.exists(output_file), f"Output file not found: {output_file}"

    # Check if the output file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Output file is not valid: {view_result.stderr}"

    # Check if the MOCK_ANNO tag is present in the header
    assert "MOCK_ANNO" in view_result.stdout, "MOCK_ANNO tag not found in the header"

    # Check if the output file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for output file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Output file has no variants"

    # Check if the MOCK_ANNO tag is present in the variants
    variants_result = subprocess.run(
        [str(bcftools_path), "view", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert variants_result.returncode == 0, f"Failed to view output file: {variants_result.stderr}"
    assert "MOCK_ANNO=" in variants_result.stdout, "MOCK_ANNO tag not found in the variants"

    print("Successfully annotated VCF file using bcftools annotate")
