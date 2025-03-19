from pathlib import Path
import subprocess
from datetime import datetime
from .base import VEPDatabase
from ..utils.validation import compute_md5


class DatabaseInitializer(VEPDatabase):
    """Handles database initialization"""
    def __init__(self, name: str, input_file: Path, fasta_ref: Path, output_dir: Path, threads: int):
        self.output_dir = Path(output_dir)
        super().__init__(self.output_dir / name)
        self.input_file = Path(input_file)
        self.fasta_ref = Path(fasta_ref)
        self.threads = max(int(threads), 1)

    def initialize(self) -> None:
        """Initialize new VEP database"""
        if not self.input_file.exists():
            raise FileNotFoundError("Input BCF file does not exist.")
        if self.variants_bcf.exists():
            raise FileExistsError("Output database already exists.")

        self.ensure_indexed(self.input_file)
        self._validate_inputs()
        self._create_database()

    def _validate_inputs(self) -> None:
        """Validate the input BCF file format"""
        is_valid, error = self.validate_bcf_header(self.input_file, norm=False)
        if not is_valid:
            raise ValueError(f"Invalid input BCF: {error}")

    def _create_database(self) -> None:
        """Create and initialize the database"""
        self.blueprint_dir.mkdir(parents=True, exist_ok=True)
        input_md5 = compute_md5(self.input_file)

        try:
            cmd = [
                "bcftools", "view", "-Ou", "--threads", str(self.threads), str(self.input_file),
                "|", "bcftools", "annotate", "-x", "INFO", "--threads", str(self.threads),
                "|", "bcftools", "norm", "-m-", "-f", str(self.fasta_ref), "-c", "x",
                "--threads", str(self.threads), "--rm-dup", "all",
                "-Ob", "--write-index", "-o", str(self.variants_bcf)
            ]

            start_time = datetime.now()
            subprocess.run(" ".join(cmd), shell=True, check=True)
            duration = datetime.now() - start_time

            # Log initialization details
            self.info_file.parent.mkdir(parents=True, exist_ok=True)
            self.info_file.write_text("")  # Create empty info file
            self.log_message(f"Initialized VEP database: {self.db_path.name}")
            self.log_message(f"Input file: {self.input_file}")
            self.log_message(f"Input file MD5: {input_md5}")
            self.log_message(f"Command executed: {' '.join(cmd)}")
            self.log_message(f"Processing completed in {duration.total_seconds():.2f} seconds")

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Error during bcftools operation: {e}")