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
        #print(f"RRRRRRRRRRRRRRRunning stash-annotate with commands: {cmd}")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result
    except Exception as e:
        print(f"Error running stash-annotate: {e}\nRuuning commands: {cmd}")
        raise e
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


def test_sample_file_validity(test_output_dir):
    """Test that the sample BCF file is valid."""
    print("\n=== Testing sample file validity ===")

    # Get bcftools path for verification
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Check if the sample BCF file exists
    assert TEST_SAMPLE.exists(), f"Sample BCF file not found: {TEST_SAMPLE}"
    print(f"Sample file exists: {TEST_SAMPLE}")

    # Check if the sample BCF file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", str(TEST_SAMPLE)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Sample BCF file is not valid: {view_result.stderr}"
    print("Sample file has valid header")

    # Check if the sample BCF file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", str(TEST_SAMPLE)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for sample BCF file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Sample BCF file has no variants"

    # Extract the number of records
    num_records = 0
    for line in stats_result.stdout.splitlines():
        if "number of records:" in line:
            num_records = int(line.split(":")[-1].strip())
            break

    print(f"Sample file has {num_records} variants")
    print("Successfully verified sample file validity")


def test_annotate_workflow(test_output_dir):
    """Test the annotate workflow using stash-annotate and annotate commands."""
    print("\n=== Testing annotate workflow ===")

    # Step 1: Run stash-init
    print("Running stash-init...")
    init_result = run_stash_init(TEST_VCF, test_output_dir, force=True)
    assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"

    # Print information about the workflow directory and files
    workflow_dir = Path(test_output_dir) / "workflow"
    print(f"Workflow directory exists: {workflow_dir.exists()}")
    if workflow_dir.exists():
        print(f"Workflow directory contents: {os.listdir(workflow_dir)}")

    # Ensure workflow directory has required files
    workflow_files = {
        'main.nf': workflow_dir / 'main.nf',
        'modules/annotate.nf': workflow_dir / 'modules' / 'annotate.nf',
        'init.yaml': workflow_dir / 'init.yaml'
    }

    # Verify all required files exist
    for name, path in workflow_files.items():
        assert path.exists(), f"Required workflow file missing: {name}"
        print(f"Found required file: {name}")

    # Step 2: Run stash-annotate with correct pathing
    print("Running stash-annotate...")
    annotate_name = "test_annotation"
    annotate_result = run_stash_annotate(test_output_dir, annotate_name, force=True)
    if annotate_result.returncode != 0:
        print(f"Command output: {annotate_result.stdout}")
        # print(f"Command error: {annotate_result.stderr}")
        print(f"Blueprint directory contents: {os.listdir(Path(test_output_dir) / 'blueprint')}")
    assert annotate_result.returncode == 0, f"stash-annotate failed: {annotate_result.stderr}"

    # Get bcftools path for verification
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Step 3: Verify the annotation directory was created
    stash_dir = os.path.join(test_output_dir, "stash")
    annotation_dir = os.path.join(stash_dir, annotate_name)
    assert os.path.exists(annotation_dir), f"Annotation directory not found: {annotation_dir}"
    print(f"Annotation directory created: {annotation_dir}")

    # Step 4: Create output directory for annotate
    output_dir = os.path.join(test_output_dir, "annotate_output")

    # Step 5: Run annotate
    print("Running annotate...")
    annotate_result = run_annotate(annotation_dir, TEST_SAMPLE, output_dir, force=True)
    if annotate_result.returncode != 0:
        print(f"Command output: {annotate_result.stdout}")
        print(f"Command error: {annotate_result.stderr}")
        print(f"Working directory contents: {os.listdir(test_output_dir)}")
        print(f"Workflow directory contents: {os.listdir(workflow_dir)}")
    assert annotate_result.returncode == 0, f"annotate failed: {annotate_result.stderr}"

    # Step 6: Verify the output directory exists
    assert os.path.exists(output_dir), f"Output directory not found: {output_dir}"
    print(f"Output directory created: {output_dir}")

    # Step 7: Verify the output file exists
    output_file = os.path.join(output_dir, os.path.splitext(os.path.basename(str(TEST_SAMPLE)))[0] + "_vst.bcf")
    if not os.path.exists(output_file):
        print(
            f"Files in output_dir ({output_dir}): {os.listdir(output_dir) if os.path.exists(output_dir) else 'Directory does not exist'}")
        raise FileNotFoundError(f"Output file not found: {output_file}")
    print(f"Output file created: {output_file}")

    # Step 8: Verify the output file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Output file is not valid: {view_result.stderr}"
    print("Output file has valid header")

    # Step 9: Verify the MOCK_ANNO tag is present in the header
    assert "MOCK_ANNO" in view_result.stdout, "MOCK_ANNO tag not found in the header"
    print("MOCK_ANNO tag found in header")

    # Step 10: Verify the output file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for output file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Output file has no variants"

    # Extract the number of records
    num_records = 0
    for line in stats_result.stdout.splitlines():
        if "number of records:" in line:
            num_records = int(line.split(":")[-1].strip())
            break

    print(f"Output file has {num_records} variants")

    # Step 11: Verify the MOCK_ANNO tag is present in the variants
    variants_result = subprocess.run(
        [str(bcftools_path), "view", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert variants_result.returncode == 0, f"Failed to view output file: {variants_result.stderr}"
    if "MOCK_ANNO=" not in variants_result.stdout:
        # Print all INFO tags present in the output
        import re
        tags = set()
        for line in variants_result.stdout.splitlines():
            if not line.startswith("#"):
                fields = line.split("\t")
                if len(fields) > 7:
                    info_field = fields[7]
                    for entry in info_field.split(";"):
                        tag = entry.split("=")[0]
                        tags.add(tag)
        print(f"Existing INFO tags in variants: {sorted(tags)}")
        assert False, "MOCK_ANNO tag not found in the variants"
    print("MOCK_ANNO tag found in variants")

    print("Successfully tested annotate workflow")


def test_annotate_with_add(test_output_dir):
    """Test the annotate workflow with stash-add using stash-annotate and annotate commands."""
    print("\n=== Testing annotate workflow with stash-add ===")

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

    # Step 5: Create output directory for annotate
    output_dir = os.path.join(test_output_dir, "annotate_output")

    # Step 6: Run annotate
    print("Running annotate...")
    annotate_result = run_annotate(annotation_dir, TEST_SAMPLE, output_dir, force=True)
    if annotate_result.returncode != 0:
        print(f"Command output: {annotate_result.stdout}")
        print(f"Command error: {annotate_result.stderr}")
        print(f"Working directory contents: {os.listdir(test_output_dir)}")
        print(f"Workflow directory contents: {os.listdir(workflow_dir)}")
    assert annotate_result.returncode == 0, f"annotate failed: {annotate_result.stderr}"

    # Step 7: Verify the output directory exists
    assert os.path.exists(output_dir), f"Output directory not found: {output_dir}"
    print(f"Output directory created: {output_dir}")

    # Step 8: Verify the output file exists
    output_file = os.path.join(output_dir, os.path.splitext(os.path.basename(str(TEST_SAMPLE)))[0] + "_vst.bcf")
    if not os.path.exists(output_file):
        print(
            f"Files in output_dir ({output_dir}): {os.listdir(output_dir) if os.path.exists(output_dir) else 'Directory does not exist'}")
        raise FileNotFoundError(f"Output file not found: {output_file}")
    print(f"Output file created: {output_file}")

    # Step 9: Verify the output file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Output file is not valid: {view_result.stderr}"
    print("Output file has valid header")

    # Step 10: Verify the MOCK_ANNO tag is present in the header
    assert "MOCK_ANNO" in view_result.stdout, "MOCK_ANNO tag not found in the header"
    print("MOCK_ANNO tag found in header")

    # Step 11: Verify the output file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for output file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Output file has no variants"

    # Extract the number of records
    num_records = 0
    for line in stats_result.stdout.splitlines():
        if "number of records:" in line:
            num_records = int(line.split(":")[-1].strip())
            break

    print(f"Output file has {num_records} variants")

    # Step 12: Verify the MOCK_ANNO tag is present in the variants
    variants_result = subprocess.run(
        [str(bcftools_path), "view", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert variants_result.returncode == 0, f"Failed to view output file: {variants_result.stderr}"
    assert "MOCK_ANNO=" in variants_result.stdout, "MOCK_ANNO tag not found in the variants"
    print("MOCK_ANNO tag found in variants")

    print("Successfully tested annotate workflow with stash-add")


def test_full_annotation_workflow(test_output_dir):
    """Test the full annotation workflow from stash-init to annotate."""
    print("\n=== Testing full annotation workflow ===")

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

    # Step 5: Create output directory for annotate
    output_dir = os.path.join(test_output_dir, "full_workflow_output")

    # Step 6: Run annotate
    print("Running annotate...")
    annotate_result = run_annotate(annotation_dir, TEST_SAMPLE, output_dir, force=True)
    if annotate_result.returncode != 0:
        print(f"Command output: {annotate_result.stdout}")
        print(f"Command error: {annotate_result.stderr}")
        print(f"Working directory contents: {os.listdir(test_output_dir)}")
        print(f"Workflow directory contents: {os.listdir(workflow_dir)}")
    assert annotate_result.returncode == 0, f"annotate failed: {annotate_result.stderr}"

    # Step 7: Verify the output directory exists
    assert os.path.exists(output_dir), f"Output directory not found: {output_dir}"
    print(f"Output directory created: {output_dir}")

    # Step 8: Verify the output file exists
    output_file = os.path.join(output_dir, os.path.splitext(os.path.basename(str(TEST_SAMPLE)))[0] + "_vst.bcf")
    if not os.path.exists(output_file):
        print(
            f"Files in output_dir ({output_dir}): {os.listdir(output_dir) if os.path.exists(output_dir) else 'Directory does not exist'}")
        raise FileNotFoundError(f"Output file not found: {output_file}")
    print(f"Output file created: {output_file}")

    # Step 9: Verify the output file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Output file is not valid: {view_result.stderr}"
    print("Output file has valid header")

    # Step 10: Verify the MOCK_ANNO tag is present in the header
    assert "MOCK_ANNO" in view_result.stdout, "MOCK_ANNO tag not found in the header"
    print("MOCK_ANNO tag found in header")

    # Step 11: Verify the output file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for output file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Output file has no variants"

    # Extract the number of records
    num_records = 0
    for line in stats_result.stdout.splitlines():
        if "number of records:" in line:
            num_records = int(line.split(":")[-1].strip())
            break

    print(f"Output file has {num_records} variants")

    # Step 12: Verify the MOCK_ANNO tag is present in the variants
    variants_result = subprocess.run(
        [str(bcftools_path), "view", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert variants_result.returncode == 0, f"Failed to view output file: {variants_result.stderr}"
    assert "MOCK_ANNO=" in variants_result.stdout, "MOCK_ANNO tag not found in the variants"
    print("MOCK_ANNO tag found in variants")

    print("Successfully tested full annotation workflow")
