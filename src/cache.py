import os
import sys
import subprocess
import argparse
from pathlib import Path
from datetime import datetime
import shutil

def log_message(log_file, message, level="INFO"):
    """Log a message with timestamp and level"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted_message = f"[{timestamp}] {level}: {message}"
    with open(log_file, "a") as log:
        log.write(formatted_message + "\n")

def check_duplicate_md5(info_file, new_md5):
    """Check if a file with the same MD5 was already added"""
    try:
        with open(info_file, "r") as f:
            for line in f:
                if "Input file MD5:" in line and new_md5 in line:
                    return True
    except FileNotFoundError:
        return False
    return False

def validate_bcf_header(bcf_path, norm:bool = True):
    """
    Validate BCF header for required normalization command and contig format.
    Returns tuple (is_valid, error_message).
    """
    try:
        header = subprocess.run(
            ["bcftools", "view", "-h", bcf_path],
            check=True,
            capture_output=True,
            text=True
        ).stdout

        if norm:
            # Check normalization command
            norm_lines = [line for line in header.splitlines()
                         if line.startswith("##bcftools_normCommand")]

            if not norm_lines:
                return False, "Missing bcftools_normCommand in header"

            norm_cmd = norm_lines[0]
            required_options = ["norm", "-c x", "-m-"]
            missing_options = [opt for opt in required_options if opt not in norm_cmd]

            if missing_options:
                return False, f"Missing required normalization options: {', '.join(missing_options)}"

        # Check contig format
        contig_lines = [line for line in header.splitlines()
                       if line.startswith("##contig=")]

        if not contig_lines:
            return False, "No contig lines found in header"

        invalid_contigs = [line for line in contig_lines
                          if not line.startswith("##contig=<ID=chr")]

        if invalid_contigs:
            example = invalid_contigs[0] if invalid_contigs else ""
            return False, f"Invalid contig format (should start with 'chr'): {example}"

        return True, None

    except subprocess.CalledProcessError as e:
        return False, f"Error reading BCF header: {e}"

def get_bcf_stats(bcf_path):
    """Get statistics from BCF file using bcftools stats"""
    try:
        result = subprocess.run(
            ["bcftools", "stats", bcf_path],
            capture_output=True,
            text=True,
            check=True
        )
        stats = {}
        for line in result.stdout.splitlines():
            if line.startswith("SN"):
                parts = line.split("\t")
                if len(parts) >= 4:
                    key = parts[2].strip(":")
                    value = parts[3]
                    stats[key] = value
        return stats
    except subprocess.CalledProcessError as e:
        return {"error": f"Failed to get statistics: {e}"}


def compute_md5(file_path):
    try:
        result = subprocess.run(
            ["md5sum", file_path],
            check=True,
            capture_output=True,
            text=True
        )
        return result.stdout.split()[0]  # The MD5 hash is the first word in the output
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error computing MD5 checksum: {e}")


def check_bcftools_installed():
    try:
        subprocess.run(["bcftools", "--version"], check=True, capture_output=True)
    except FileNotFoundError:
        sys.exit("Error: bcftools is not installed or not in PATH.")


def log_script_command(info_file):
    """Log the exact command used to execute the script"""
    command = f"python3 {sys.argv[0]} {' '.join(sys.argv[1:])}"
    log_message(info_file, f"Script command: {command}")

def store_workflow_dag(run_dir:Path, cmd:list):
    """Generate and store workflow DAG visualization"""

    try:
        del cmd[cmd.index("-with-trace")]
    except ValueError:
        pass  # Element not in list

    try:
        # Generate DAG in dot format
        dag_cmd = cmd + [ "-preview", "-with-dag", str(run_dir / "flowchart.html") ]

        subprocess.run(dag_cmd, check=True, capture_output=True)


    except subprocess.CalledProcessError as e:
        log_message(run_dir / "annotation.info", f"Warning: Failed to generate workflow DAG: {e}", level="WARN")


def init_mode(args):

    # Add validation at start, will be normalized anyway
    is_valid, error = validate_bcf_header(args.i, norm=False)
    if not is_valid:
        sys.exit(f"Invalid input BCF: {error}")

    output_dir = Path(args.name)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_bcf = output_dir / "vep_db.bcf"
    info_file = output_dir / "vep_db.bcf.info"

    if not Path(args.i).exists():
        sys.exit("Error: Input BCF file does not exist.")
    if output_bcf.exists():
        sys.exit("Error: Output database already exists.")

    threads = str(max(int(args.t), 1))
    input_md5 = compute_md5(args.i)

    try:
        if args.sort_mem:
            cmd = [
                "bcftools", "annotate", "-x", "INFO", "--threads", threads, args.i,
                "|", "bcftools", "norm", "-m-", "-f", args.fasta, "-c", "x", "--threads", threads, "--rm-dup", "all",
                "|", "bcftools", "sort", "--threads", threads, "-m", args.sort_mem, "-Ob", "--write-index", "-o",
                str(output_bcf)
            ]
        else:
            cmd = [
                "bcftools", "annotate", "-x", "INFO", "--threads", threads, args.i,
                "|", "bcftools", "norm", "-m-", "-f", args.fasta, "-c", "x", "--threads", threads, "--rm-dup", "all",
                "-Ob", "--write-index", "-o", str(output_bcf)
            ]

        start_time = datetime.now()
        subprocess.run(" ".join(cmd), shell=True, check=True)
        duration = datetime.now() - start_time

        # Only log after successful creation of BCF file
        info_file.write_text("")  # Create empty info file
        log_script_command(info_file)
        log_message(info_file, f"Initialized VEP database: {args.name}")
        log_message(info_file, f"Input file: {args.i}")
        log_message(info_file, f"Input file MD5: {input_md5}")
        log_message(info_file, f"Command executed: {' '.join(cmd)}")
        log_message(info_file, f"Processing completed in {duration.total_seconds():.2f} seconds")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error during bcftools operation: {e}")

def replace_db(db_bcf, temp_output_bcf):
    # Add validation for both files
    is_valid, error = validate_bcf_header(db_bcf)
    if not is_valid:
        sys.exit(f"Invalid input BCF: {error}")
    is_valid, error = validate_bcf_header(temp_output_bcf)
    if not is_valid:
        sys.exit(f"Invalid database BCF: {error}")

    # Get statistics before replacement
    info_file = Path(db_bcf).parent / "vep_db.bcf.info"
    stats = get_bcf_stats(temp_output_bcf)
    new_md5 = compute_md5(temp_output_bcf)

    # Replace files
    os.replace(temp_output_bcf, db_bcf)
    os.replace(str(temp_output_bcf) + '.csi', str(db_bcf) + '.csi')

    # Log statistics
    log_message(info_file, "Database statistics after update:")
    log_message(info_file, f"Database MD5: {new_md5}")
    if "number of SNPs:" in stats:
        log_message(info_file, f"Number of SNPs: {stats['number of SNPs']}")
    if "number of indels:" in stats:
        log_message(info_file, f"Number of indels: {stats['number of indels']}")
    if "number of records:" in stats:
        log_message(info_file, f"Total variants: {stats['number of records']}")
    if "number of samples:" in stats:
        log_message(info_file, f"Number of samples: {stats['number of samples']}")


def add_mode(args):

    db_bcf = Path(args.db) / "vep_db.bcf"
    info_file = Path(args.db) / "vep_db.bcf.info"
    input_vcf = Path(args.vcf)
    threads = str(max(int(args.threads), 1))

    if not db_bcf.exists() or not input_vcf.exists():
        sys.exit("Error: Database or input VCF file does not exist.")

    input_md5 = compute_md5(input_vcf)
    if check_duplicate_md5(info_file, input_md5):
        sys.exit("Error: This file was already added to the database (MD5 match).")

    try:
        temp_output_bcf = db_bcf.parent / "temp_vep_db.bcf"
        start_time = datetime.now()

        if args.skip_preprocessing:

            cmd_merge = ["bcftools", "merge", "-m", "none", "-o", str(temp_output_bcf), "-Ob", "--write-index", "--threads", threads,
                         str(db_bcf), str(input_vcf)]
            subprocess.run(cmd_merge, check=True)
            duration = datetime.now() - start_time
            # we might need this afterwards?
            # bcftools norm vep_db.bcf -m- -c x --threads 10 -f /mnt/data/resources/reference/ucsc/hg19_canonical.fa.gz |
            # bcftools sort -Ob -o vep_db_s.bcf --write-index -m 10G -T tmp/

            # Replace the old database with the new one
            replace_db(db_bcf, temp_output_bcf)

            # Log only after successful operation
            log_script_command(info_file)
            log_message(info_file, f"Added new VCF: {input_vcf}")
            log_message(info_file, f"Input file MD5: {input_md5}")
            log_message(info_file, f"Merge command: {' '.join(cmd_merge)}")
            log_message(info_file, f"Merge completed in {duration.total_seconds():.2f} seconds")
        else:
            temp_dir = db_bcf.parent / "add_temp"
            temp_dir.mkdir(parents=True, exist_ok=True)
            processed_input = temp_dir / "processed_input.bcf"

            try:
                if args.sort_mem:
                    preprocess_cmd = [
                        "bcftools", "annotate", "-x", "INFO", "--threads", threads, str(input_vcf),
                        "|", "bcftools", "norm", "-m-", "-f", args.fasta, "-c", "x", "--threads", threads, "--rm-dup", "all",
                        "|", "bcftools", "sort", "-m", args.sort_mem, "-Ob", "--write-index",
                        "-o", str(processed_input), "-T", str(temp_dir), "--threads", threads
                    ]
                else:
                    preprocess_cmd = [
                        "bcftools", "annotate", "-x", "INFO", "--threads", threads, str(input_vcf),
                        "|", "bcftools", "norm", "-m-", "-f", args.fasta, "-c", "x", "--threads", threads,
                        "-Ob", "--write-index", "-o", str(processed_input), "--rm-dup", "all",
                    ]

                subprocess.run(" ".join(preprocess_cmd), shell=True, check=True)
                preprocess_duration = datetime.now() - start_time

                merge_start = datetime.now()
                # Merge with specified chromosome order
                cmd_merge = ["bcftools", "merge", "-m", "none",
                             "-o", str(temp_output_bcf),
                             "-Ob", "--write-index",
                             "--threads", threads,
                             str(db_bcf), str(input_vcf)]
                subprocess.run(cmd_merge, check=True)
                merge_duration = datetime.now() - merge_start

                # Replace the old database with the new one
                replace_db(db_bcf, temp_output_bcf)

                # Log only after successful operation
                log_script_command(info_file)
                log_message(info_file, f"Added new VCF: {input_vcf}")
                log_message(info_file, f"Input file MD5: {input_md5}")
                log_message(info_file, f"Preprocessing command: {' '.join(preprocess_cmd)}")
                log_message(info_file, f"Preprocessing completed in {preprocess_duration.total_seconds():.2f} seconds")
                log_message(info_file, f"Merge command: {' '.join(cmd_merge)}")
                log_message(info_file, f"Merge completed in {merge_duration.total_seconds():.2f} seconds")
            finally:
                import shutil
                shutil.rmtree(temp_dir)

        log_message(info_file, "Database update completed successfully")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error during bcftools operation: {e}")
    except Exception as e:
        raise


def create_unique_annotation_dir(db_dir, workflow_hash) -> Path:
    """Create a unique directory for this annotation run using timestamp and workflow hash"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    annotations_dir = Path(db_dir) / "annotations"
    run_dir = annotations_dir / f"{timestamp}_{workflow_hash[:8]}"  # Use first 8 chars of hash
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def get_workflow_hash(workflow_dir):
    """Get combined hash of workflow files"""
    workflow_files = ['main.nf', 'nextflow.config']
    combined_content = b''

    for file in workflow_files:
        path = Path(workflow_dir) / file
        if not path.exists():
            raise FileNotFoundError(f"Required workflow file not found: {path}")
        combined_content += path.read_bytes()

    import hashlib
    return hashlib.md5(combined_content).hexdigest()


def store_workflow_files(workflow_dir, target_dir):
    """Store workflow files with their hash for version control"""
    workflow_files = {
        'main.nf': workflow_dir / 'main.nf',
        'nextflow.config': workflow_dir / 'nextflow.config'
    }

    stored_files = {}
    for name, path in workflow_files.items():
        if not path.exists():
            raise FileNotFoundError(f"Required workflow file not found: {path}")

        # Copy file
        target_path = target_dir / name
        with open(path, 'rb') as src, open(target_path, 'wb') as dst:
            dst.write(src.read())

        # Store hash
        stored_files[name] = compute_md5(target_path)

    return stored_files


def find_latest_annotation_run(db_dir):
    """Find the latest annotation run directory"""
    annotations_dir = Path(db_dir) / "annotations"
    if not annotations_dir.exists():
        return None

    runs = [d for d in annotations_dir.glob("*_*") if d.is_dir()]
    if not runs:
        return None

    return max(runs, key=lambda x: x.stat().st_mtime)


def validate_workflow_version(annotations_dir, workflow_hashes):
    """Validate that workflow version matches previous runs"""
    latest_run = find_latest_annotation_run(annotations_dir)
    if not latest_run:
        return True, None  # No previous runs, so valid

    latest_info = latest_run / "annotation.info"
    try:
        with open(latest_info, 'r') as f:
            content = f.read()
            for name, hash_value in workflow_hashes.items():
                expected_line = f"Workflow {name} MD5: {hash_value}"
                if expected_line not in content:
                    return False, f"Workflow file {name} has changed since last annotation run"
        return True, None
    except FileNotFoundError:
        return False, "Previous annotation info file not found"


def annotate_mode(args):
    db_bcf = Path(args.db) / "vep_db.bcf"
    db_info = Path(args.db) / "vep_db.bcf.info"

    if not db_bcf.exists():
        sys.exit("Error: Database BCF file does not exist.")

    # Get workflow hash and create directory
    workflow_hash = get_workflow_hash(args.workflow)
    run_dir = create_unique_annotation_dir(args.db, workflow_hash)

    # Create run-specific info file
    run_info = run_dir / "annotation.info"

    # Copy original info file
    shutil.copy2(db_info, run_dir / "vep_db.bcf.info.snapshot")

    # Store workflow files and get their hashes
    workflow_files = store_workflow_files(Path(args.workflow), run_dir)

    # For subsequent runs, validate workflow version unless force flag is used
    previous_runs = list(Path(args.db).glob("annotations/*"))
    if previous_runs and not args.force:
        is_valid, error = validate_workflow_version(args.db, workflow_files)
        if not is_valid:
            sys.exit(f"Workflow validation failed: {error}\nUse --force to override this check.")


    try:
        # Run nextflow workflow in database mode first
        cmd = [
            "nextflow", "run", str(run_dir / "main.nf"),
            "-with-trace",
            "--input", str(db_bcf),
            "--output", str(run_dir),
            "--db_mode", "true"  # Enable database mode
        ]
        if args.params:
            cmd.extend(["--params-file", args.params])

        start_time = datetime.now()
        subprocess.run(cmd, check=True)

        # Store workflow DAG
        store_workflow_dag(run_dir, cmd)

        duration = datetime.now() - start_time

        # Log annotation details
        log_message(run_info, f"Annotated database: {db_bcf}")
        log_message(run_info, f"Workflow hash: {workflow_hash}")
        log_message(run_info, f"Annotation directory: {run_dir}")
        log_message(run_info, f"Command executed: {' '.join(cmd)}")
        for name, hash_value in workflow_files.items():
            log_message(run_info, f"Workflow {name} MD5: {hash_value}")
        log_message(run_info, f"Processing completed in {duration.total_seconds():.2f} seconds")

        # Store tool versions in run directory
        tool_versions = Path(run_dir) / "tool_version.log"
        if tool_versions.exists():
            with open(tool_versions, 'r') as f:
                for line in f:
                    log_message(run_info, f"Tool version: {line.strip()}")

    except subprocess.CalledProcessError as e:
        sys.exit(f"Error during workflow execution: {e}")


def main():
    check_bcftools_installed()

    parser = argparse.ArgumentParser(
        description="Manage VEP database with BCFtools."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # init command
    init_parser = subparsers.add_parser("init", help="Initialize VEP database")
    init_parser.add_argument("-n", "--name", dest="name", required=True, help="Name of the new database directory")
    init_parser.add_argument("-i", "--vcf", dest="i", help="CSI-indexed BCF file")
    init_parser.add_argument("-p", "--parent-dir", dest="p", default=".",
                            help="Parent output directory (default: current directory)")
    init_parser.add_argument("-t", "--threads", dest="t", help="Use multithreading with <int> worker threads [4]", default=4)
    init_parser.add_argument("-f", "--fasta", required=True, help="Reference FASTA file for variant normalization")
    init_parser.add_argument("-s", "--sort-mem", help="Sort output with specified memory limit (e.g., '768M'). Skip sorting if not provided")

    # add command
    add_parser = subparsers.add_parser("add", help="Add new VCF to the database")
    add_parser.add_argument("-d", "--db", required=True, help="Path to the existing database directory")
    add_parser.add_argument("-i", "--vcf", help="Path to the VCF file to be added")
    add_parser.add_argument("-t", "--threads", help="Use multithreading with <int> worker threads [4]", default=4)
    add_parser.add_argument("-f", "--fasta", required=True, help="Reference FASTA file for variant normalization")
    add_parser.add_argument("-s", "--sort-mem", help="Sort output with specified memory limit (e.g., '768M'). Skip sorting if not provided")
    add_parser.add_argument("--skip-preprocessing", action="store_true", help="Skip preprocessing steps (annotate, norm, sort) and merge input directly")

    # annotate command
    annotate_parser = subparsers.add_parser("annotate", help="Run annotation workflow on database")
    annotate_parser.add_argument("-d", "--db", required=True, help="Path to the existing database directory")
    annotate_parser.add_argument("-w", "--workflow", required=True, help="Directory containing workflow files (main.nf, nextflow.config)")
    annotate_parser.add_argument("-p", "--params", help="Optional parameters file for nextflow workflow")
    annotate_parser.add_argument("-f", "--force", action="store_true", help="Force annotation even if workflow files have changed")

    args = parser.parse_args()

    if args.command == "init":
        # Create info file first so we can log the command
        output_dir = Path(args.name)
        output_dir.mkdir(parents=True, exist_ok=True)
        info_file = output_dir / "vep_db.bcf.info"
        info_file.write_text("")
        log_script_command(info_file)
        init_mode(args)
    elif args.command == "add":
        # Log command to existing info file
        info_file = Path(args.db) / "vep_db.bcf.info"
        log_script_command(info_file)
        add_mode(args)
    elif args.command == "annotate":
        info_file = Path(args.db) / "vep_db.bcf.info"
        log_script_command(info_file)
        annotate_mode(args)

if __name__ == "__main__":
    main()
