#!/usr/bin/env python3
"""
update_reference_data.py - Script to update reference data for tests
"""

import re
import sys
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import uuid

TEST_ROOT = os.path.dirname(os.path.abspath(__file__))
VCFSTASH_CMD = os.path.join(os.path.dirname(TEST_ROOT), "vcfstash.py")
TEST_DATA_DIR = os.path.join(TEST_ROOT, "data", "nodata")
TEST_CONFIG = os.path.join(TEST_ROOT, "config", "nextflow_test.config")
TEST_PARAMS = os.path.join(TEST_ROOT, "config", "user_params.yaml")
TEST_VCF = os.path.join(TEST_DATA_DIR, "crayz_db.bcf")
EXPECTED_OUTPUT_DIR = os.path.join(TEST_ROOT, "data", "expected_output")
TEST_ANNO_CONFIG = os.path.join(os.path.dirname(__file__), "config", "annotation.config")


def update_stash_add_reference_data(force=True):
    """Update the reference data for stash-add function."""
    # Define path constants

    REFERENCE_DIR = os.path.join(EXPECTED_OUTPUT_DIR, "stash_add_annotate_result")

    print(f"=== Updating reference data for stash-add ===")
    print(f"Target directory: {REFERENCE_DIR}")

    # First, ensure stash-init data exists
    init_reference_dir = os.path.join(EXPECTED_OUTPUT_DIR, "stash_init_result")
    if not os.path.exists(init_reference_dir):
        print("Failed to create stash-init data, cannot proceed with stash-add")
        return None

    # Create a temporary directory for stash-add
    temp_dir = tempfile.mkdtemp(prefix="vcfstash_ref_")

    try:
        # Copy the existing stash-init reference to our temp directory
        print(f"Copying stash-init reference data to {temp_dir}")
        for item in os.listdir(init_reference_dir):
            src = os.path.join(init_reference_dir, item)
            dst = os.path.join(temp_dir, item)
            if os.path.isdir(src):
                shutil.copytree(src, dst)
            else:
                shutil.copy2(src, dst)

        # Define the second VCF file
        test_vcf2 = os.path.join(TEST_DATA_DIR, "crayz_db2.bcf")

        # Verify the second VCF file exists
        if not os.path.exists(test_vcf2):
            print(f"Error: Second VCF file not found at {test_vcf2}")
            return None

        # Run stash-add
        add_cmd = [
            VCFSTASH_CMD,
            "stash-add",
            "--db", temp_dir,
            "-i", test_vcf2
        ]

        print(f"Running command: {' '.join(add_cmd)}")

        add_result = subprocess.run(
            add_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if add_result.returncode != 0:
            print(f"stash-add failed with exit code {add_result.returncode}")
            print(f"stdout: {add_result.stdout}")
            print(f"stderr: {add_result.stderr}")
            return None

        print(f"stash-add completed successfully")

        # Create or clear the reference directory
        if os.path.exists(REFERENCE_DIR):
            shutil.rmtree(REFERENCE_DIR)
        os.makedirs(REFERENCE_DIR, exist_ok=True)

        # Copy all files from temp_dir to reference_dir
        print(f"Copying files from {temp_dir} to {REFERENCE_DIR}")
        for item in os.listdir(temp_dir):
            source = os.path.join(temp_dir, item)
            dest = os.path.join(REFERENCE_DIR, item)

            if os.path.isdir(source):
                shutil.copytree(source, dest)
            else:
                shutil.copy2(source, dest)

        print("Reference data for stash-add updated successfully")
        print(f"Output saved to: {REFERENCE_DIR}")
        return REFERENCE_DIR

    except Exception as e:
        print(f"Error during reference data update: {str(e)}")
        import traceback
        traceback.print_exc()
        return None
    finally:
        # Clean up the temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)


def normalize_bcf_timestamps(bcf_file):
    """Normalize timestamps in BCF file to make tests more stable."""
    if not os.path.exists(bcf_file):
        print(f"Warning: BCF file not found at {bcf_file}")
        return

    print(f"Normalizing timestamps in {bcf_file}")

    # Get bcftools path
    bcftools_path = os.path.join(os.environ.get('VCFSTASH_ROOT', ''), 'tools', 'bcftools')
    if not os.path.exists(bcftools_path):
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Create a temporary VCF file
    temp_vcf = bcf_file + ".temp.vcf"

    # Convert BCF to VCF
    result = subprocess.run(
        [bcftools_path, "view", bcf_file, "-o", temp_vcf],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    if result.returncode != 0:
        print(f"Error converting BCF to VCF: {result.stderr}")
        if os.path.exists(temp_vcf):
            os.remove(temp_vcf)
        return

    # Read and modify the VCF content
    modified_lines = []
    with open(temp_vcf, 'r') as f:
        for line in f:
            # Replace timestamp patterns in header lines
            if line.startswith('##'):
                line = re.sub(r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d+\s+\d+:\d+:\d+\s+\d+\b',
                              'Jan 1 00:00:00 2023', line)
                # Also handle fileDate format
                if "##fileDate=" in line:
                    line = "##fileDate=20230101\n"
            modified_lines.append(line)

    # Write modified content back to the temporary file
    with open(temp_vcf, 'w') as f:
        f.writelines(modified_lines)

    # Convert back to BCF
    result = subprocess.run(
        [bcftools_path, "view", "-O", "b", "-o", bcf_file, temp_vcf],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    if result.returncode != 0:
        print(f"Error converting normalized VCF back to BCF: {result.stderr}")

    # Clean up
    if os.path.exists(temp_vcf):
        os.remove(temp_vcf)


def normalize_text_file_timestamps(file_path):
    """Normalize timestamps in text files to make tests more stable."""
    if not os.path.exists(file_path):
        print(f"Warning: File not found at {file_path}")
        return

    print(f"Normalizing timestamps in {file_path}")

    # Read the file content
    with open(file_path, 'r') as f:
        content = f.read()

    # Replace timestamp patterns
    # ISO format: 2023-04-02T18:14:50
    normalized_content = re.sub(
        r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}',
        '2023-01-01T00:00:00',
        content
    )

    # Date formats like "Apr 2 18:14:50 2025"
    normalized_content = re.sub(
        r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d+\s+\d+:\d+:\d+\s+\d+\b',
        'Jan 1 00:00:00 2023',
        normalized_content
    )

    # Write normalized content back
    with open(file_path, 'w') as f:
        f.write(normalized_content)


# Add to update_reference.py
def update_annotation_reference_data(force: bool = False) -> bool:
    """Update the reference data for testing the annotation functionality.

    Creates a reference database with annotations that can be used in tests.
    The reference data includes normalized timestamps and standardized configuration.

    Args:
        force: If True, overwrite existing reference data. Defaults to False.

    Returns:
        bool: True if reference data was successfully updated, False otherwise.
    """ 

    # Reference annotation name
    annotation_name = "test_annotation"

    # Define the base output directory
    annotation_ref_dir = os.path.join(EXPECTED_OUTPUT_DIR, "stash_add_annotate_result",
                                      "stash", annotation_name)

    print(f"Using annotation reference directory: {annotation_ref_dir}")

    # Check if we already have reference data and not forcing an update
    if os.path.exists(annotation_ref_dir) and not force:
        print(f"Annotation reference data exists at {annotation_ref_dir}. Use --force to overwrite.")
        return True

    # Create a parent temporary directory
    temp_parent_dir = tempfile.mkdtemp(prefix="vcfstash_anno_parent_")
    try:
        # Create a unique directory name for the stash database that doesn't exist yet
        unique_id = str(uuid.uuid4())
        tmp_dir = os.path.join(temp_parent_dir, f"stash_db_{unique_id}")

        # Ensure the directory doesn't exist
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

        # Initialize a stash with a non-existent output directory
        init_cmd = [
            sys.executable, VCFSTASH_CMD,
            "stash-init",
            "-i", TEST_VCF,
            "-o", tmp_dir,
            "-y", TEST_PARAMS,
            "-f"  # Force flag
        ]

        print(f"Running stash-init: {' '.join(init_cmd)}")
        try:
            result = subprocess.run(init_cmd, check=True, capture_output=True, text=True)
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"stash-init failed: {e}")
            print(f"stdout: {e.stdout}")
            print(f"stderr: {e.stderr}")
            return False

        # Copy the existing annotation config
        annotation_config = os.path.join(temp_parent_dir, "annotation.config")
        shutil.copyfile(TEST_ANNO_CONFIG, annotation_config)
        print(f"Copied annotation config from {TEST_ANNO_CONFIG} to {annotation_config}")

        print(f"Created annotation config at {annotation_config}")

        # Run annotation with the config in the correct format
        annotate_cmd = [
            sys.executable, VCFSTASH_CMD,
            "stash-annotate",
            "--name", annotation_name,
            "-a", annotation_config,
            "--db", tmp_dir,
            "-f"  # Force flag
        ]

        print(f"Running stash-annotate: {' '.join(annotate_cmd)}")
        try:
            result = subprocess.run(annotate_cmd, check=True, capture_output=True, text=True, timeout=300)
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"stash-annotate failed: {e}")
            print(f"stdout: {e.stdout}")
            print(f"stderr: {e.stderr}")
            return False
        except subprocess.TimeoutExpired:
            print("stash-annotate timed out after 5 minutes")
            return False

        # Make sure the expected paths exist
        tmp_annotation_dir = os.path.join(tmp_dir, "stash", annotation_name)
        if not os.path.exists(tmp_annotation_dir):
            print(f"Expected annotation directory was not created: {tmp_annotation_dir}")
            return False

        # Create the target directory structure if it doesn't exist
        os.makedirs(os.path.dirname(annotation_ref_dir), exist_ok=True)

        # Remove the existing reference directory if it exists
        if os.path.exists(annotation_ref_dir):
            shutil.rmtree(annotation_ref_dir)

        # Copy the temporary annotation directory to the reference location
        shutil.copytree(tmp_annotation_dir, annotation_ref_dir)

        # Normalize the timestamps in BCF files
        for root, _, files in os.walk(annotation_ref_dir):
            for file in files:
                if file.endswith('.bcf'):
                    bcf_path = os.path.join(root, file)
                    normalize_bcf_timestamps(bcf_path)
                elif file.endswith('.html') or file.endswith('.txt') or file.endswith('.config'):
                    text_path = os.path.join(root, file)
                    normalize_text_file_timestamps(text_path)
                elif file == 'blueprint_snapshot.info':
                    # Handle the JSON snapshot file to remove/normalize timestamps
                    import json
                    snapshot_path = os.path.join(root, file)
                    try:
                        with open(snapshot_path, 'r') as f:
                            snapshot_data = json.load(f)

                        # Normalize or remove timestamp information
                        if 'timestamp' in snapshot_data:
                            snapshot_data['timestamp'] = "NORMALIZED_TIMESTAMP"

                        # Write back the normalized data
                        with open(snapshot_path, 'w') as f:
                            json.dump(snapshot_data, f, indent=2)
                    except Exception as e:
                        print(f"Error normalizing snapshot file: {e}")

        print(f"Successfully updated annotation reference data in {annotation_ref_dir}")
        return True

    finally:
        # Clean up the temporary parent directory
        if os.path.exists(temp_parent_dir):
            shutil.rmtree(temp_parent_dir)

def update_annotate_input(force: bool = False) -> bool:
    """Store reference files for annotate command.

    Args:
        force: Whether to overwrite existing reference files

    Returns:
        bool: True if reference files were successfully updated, False otherwise.
    """
    test_input = Path(TEST_DATA_DIR) / "sample4.bcf"
    test_cache = Path(EXPECTED_OUTPUT_DIR) / "stash_add_annotate_result/stash/test_annotation"
    annotate_result = Path(EXPECTED_OUTPUT_DIR) / "annotate_result"
    cached_dir = annotate_result / "cached"
    uncached_dir = annotate_result / "uncached"

    # Verify test input and cache exist
    if not test_input.exists():
        print(f"Test input file not found: {test_input}")
        return False

    if not test_cache.exists():
        print(f"Test cache directory not found: {test_cache}")
        return False

    if annotate_result.exists():
        if force:
            print(f"Removing existing reference directory: {annotate_result}")
            shutil.rmtree(annotate_result)
        else:
            print(f"Reference directory {annotate_result} already exists. Use --force to overwrite.")
            return True

    # Create reference directories
    print(f"Creating reference directory: {annotate_result}")
    annotate_result.mkdir(parents=True, exist_ok=True)

    # Run cached annotation
    print("Running cached annotation...")
    cached_cmd = [
        sys.executable, VCFSTASH_CMD,
        "annotate",
        "-i", str(test_input),
        "-a", str(test_cache),
        "-o", str(cached_dir),
        "-f"
    ]

    try:
        subprocess.run(cached_cmd, check=True, capture_output=True, text=True, timeout=300)
    except subprocess.CalledProcessError as e:
        print(f"Cached annotation failed: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return False
    except subprocess.TimeoutExpired:
        print("Cached annotation timed out after 5 minutes")
        return False

    # Verify cached output exists
    if not cached_dir.exists() or not list(cached_dir.glob("*.bcf")):
        print(f"Cached annotation did not produce expected output in {cached_dir}")
        return False

    # Run uncached annotation
    print("Running uncached annotation...")
    uncached_cmd = [
        sys.executable, VCFSTASH_CMD,
        "annotate",
        "-i", str(test_input),
        "-a", str(test_cache),
        "-o", str(uncached_dir),
        "--uncached",
        "-f"
    ]

    try:
        subprocess.run(uncached_cmd, check=True, capture_output=True, text=True, timeout=300)
    except subprocess.CalledProcessError as e:
        print(f"Uncached annotation failed: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return False
    except subprocess.TimeoutExpired:
        print("Uncached annotation timed out after 5 minutes")
        return False

    # Verify uncached output exists
    if not uncached_dir.exists() or not list(uncached_dir.glob("*.bcf")):
        print(f"Uncached annotation did not produce expected output in {uncached_dir}")
        return False

    print("Reference files stored successfully")
    return True

# Update the main part of the script to include the new function
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Update reference data for vcfstash tests")

    parser.add_argument('--force', action='store_true', help='Force overwrite of existing reference data')
    parser.add_argument('--all', action='store_true', help='Update all reference data')
    parser.add_argument('--init', action='store_true', help='Update stash-init reference data')
    parser.add_argument('--add', action='store_true', help='Update stash-add reference data')
    parser.add_argument('--annotate', action='store_true', help='Update stash-annotate reference data')
    parser.add_argument('--annotate-input', dest="annotate_input", action='store_true', help='Update stash-annotate reference data')

    args = parser.parse_args()

    # Update stash-init reference data
    if args.init or args.all:
        success = update_reference_data(force=args.force)
        if not success:
            print("Failed to update stash-init reference data")

    # Update stash-add reference data
    if args.add or args.all:
        success = update_stash_add_reference_data(force=args.force)
        if not success:
            print("Failed to update stash-add reference data")

    # Update stash-annotate reference data
    if args.annotate or args.all:
        success = update_annotation_reference_data(force=args.force)
        if not success:
            print("Failed to update stash-annotate reference data")

    if args.annotate_input or args.all:
        success = update_annotate_input(args.force)
        if not success:
            print("Failed to update annotate input vcf reference data")
