import subprocess
import sys
from pathlib import Path
from typing import Tuple, List, Optional

import yaml
import hashlib
from src.utils.validation import ensure_indexed, validate_bcf_header
from src.utils.logging import log_message

class VEPDatabase:
    """Base class for VEP database operations"""
    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        self.blueprint_dir = self.db_path / "blueprint"
        self.variants_bcf = self.blueprint_dir / "variants.bcf"
        self.info_file = self.blueprint_dir / "variants.info"

    def ensure_indexed(self, file_path: Path) -> None:
        """Validate that BCF/VCF file has an index"""
        ensure_indexed(file_path)

    def validate_bcf_header(self, bcf_path: Path, norm: bool = True) -> Tuple[bool, str]:
        """Validate BCF header format"""
        return validate_bcf_header(bcf_path, norm)

    def log_message(self, message: str, level: str = "INFO") -> None:
        """Log a message to the database info file"""
        log_message(self.info_file, message, level)

class NextflowWorkflow:
    """Base class for Nextflow workflow operations"""
    def __init__(self, workflow_dir: Path, filename: str = "main", params_file: Optional[Path] = None):
        self.workflow_dir = Path(workflow_dir)
        self.workflow_path = self.workflow_dir / f"{filename}.nf"
        if not self.workflow_path.exists():
            raise FileNotFoundError(f"Workflow file not found: {self.workflow_path}")
        self.workflow_config_path = self.workflow_dir / "nextflow.config"
        if not self.workflow_config_path.exists():
            raise FileNotFoundError(f"Nextflow config file not found: {self.workflow_config_path}")
        self.params_file = Path(params_file)
        if self.params_file and not self.params_file.exists():
            raise FileNotFoundError(f"Parameters file not found: {self.params_file}")
        self.workflow_hash = self.get_workflow_hash(self.workflow_dir)
        self.config = self.load_nextflow_config()
        self.run_dir = self.workflow_dir / "run"
        self.run_dir.mkdir(parents=True, exist_ok=True)


    def load_nextflow_config(self, test_mode=False):
        """Load Nextflow configuration from YAML files."""

        # Load and return the configuration
        with open(self.params_file) as f:
            return yaml.safe_load(f)

    def get_workflow_hash(self, workflow_dir: Path) -> str:
        """Get combined hash of workflow files"""
        workflow_files = ['main.nf', 'nextflow.config']
        combined_content = b''

        for file in workflow_files:
            path = Path(workflow_dir) / file
            if not path.exists():
                raise FileNotFoundError(f"Required workflow file not found: {path}")
            combined_content += path.read_bytes()
        return hashlib.md5(combined_content).hexdigest()

    def store_workflow_dag(self, run_dir: Path, cmd: List[str]) -> None:
        """Generate and store workflow DAG visualization"""

        try:
            del cmd[cmd.index("-with-trace")]
        except ValueError:
            pass  # Element not in list

        try:
            # Generate DAG in dot format
            dag_cmd = cmd + ["-preview", "-with-dag", str(run_dir / "flowchart.html")]

            subprocess.run(dag_cmd, check=True, capture_output=True)


        except subprocess.CalledProcessError as e:
            log_message(run_dir / "annotation.info", f"Warning: Failed to generate workflow DAG: {e}", level="WARN")



    def run(self, input_file: Path, output_dir: Path, db_mode: bool = False, **kwargs) -> subprocess.CompletedProcess:
        """Run Nextflow pipeline with the appropriate configuration.

            Args:
                input_file: Path to input BCF/VCF file
                output_dir: Path to output directory
                db_mode: Whether to run in database mode
                **kwargs: Additional arguments to pass to the workflow
            """

        if not self.workflow_path.is_file():
            raise FileNotFoundError(f"Required workflow file not found: {self.workflow_path}")

        cmd = [
            "nextflow", "run", str(self.workflow_path),
            "-with-trace",
            "--input", str(input_file),
            "--output", str(output_dir)
        ]

        if db_mode:
            cmd.extend(["--db_mode", "true"])

        if self.params_file:
            cmd.extend(["-params-file", str(self.params_file)])

        # Add any additional arguments, filtering out empty ones
        for key, value in kwargs.items():
            if key == "nextflow_args" and value:
                # Filter out empty or None values
                args = [arg for arg in value if arg and arg != '[]']
                if args:
                    cmd.extend(args)
            elif value:  # Only add if value is not empty
                cmd.extend([f"--{key}", str(value)])

        # Print command
        print("Executing Nextflow command:", file=sys.stderr)
        print(" ".join(str(x) for x in cmd), file=sys.stderr)
        print("\n", file=sys.stderr)

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(result.stdout, file=sys.stderr)
            self.store_workflow_dag(output_dir, cmd)
            return result
        except subprocess.CalledProcessError as e:
            print("Nextflow execution failed:", file=sys.stderr)
            print("Exit code:", e.returncode, file=sys.stderr)
            print("\nSTDOUT:", file=sys.stderr)
            print(e.stdout, file=sys.stderr)
            print("\nSTDERR:", file=sys.stderr)
            print(e.stderr, file=sys.stderr)
            raise

