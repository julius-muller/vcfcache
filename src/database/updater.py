from pathlib import Path
import subprocess
import os
from datetime import datetime
from .base import VEPDatabase
from ..utils.validation import compute_md5, get_bcf_stats, check_duplicate_md5

class DatabaseUpdater(VEPDatabase):
    """Handles adding new variants to database"""
    def __init__(self, db_path: Path, input_file: Path, fasta_ref: Path, threads: int):
        super().__init__(db_path)
        self.input_file = Path(input_file)
        self.fasta_ref = Path(fasta_ref)
        self.threads = max(int(threads), 1)

    def add(self) -> None:
        """Add new variants to existing database"""
        self._validate_inputs()
        self._merge_variants()

    def _validate_inputs(self) -> None:
        """Validate input files and check for duplicates"""
        if not self.variants_bcf.exists():
            raise FileNotFoundError("Database BCF file does not exist")
        if not self.input_file.exists():
            raise FileNotFoundError("Input VCF/BCF file does not exist")

        self.ensure_indexed(self.variants_bcf)
        self.ensure_indexed(self.input_file)

        input_md5 = compute_md5(self.input_file)
        if check_duplicate_md5(self.info_file, input_md5):
            raise ValueError("This file was already added to the database (MD5 match)")

    def _merge_variants(self) -> None:
        """Merge new variants into the database"""
        temp_merged = Path(str(self.variants_bcf) + ".tmp.bcf")

        try:
            start_time = datetime.now()

            # Merge variants using bcftools
            subprocess.run([
                "bcftools", "concat",
                "--allow-overlaps",
                "--rm-dup", "all",
                "-Ob",
                "--write-index",
                "-o", str(temp_merged),
                "--threads", str(self.threads),
                str(self.variants_bcf),
                str(self.input_file)
            ], check=True)

            duration = datetime.now() - start_time
            stats = get_bcf_stats(temp_merged)

            # Replace old database with merged file
            os.replace(temp_merged, self.variants_bcf)
            os.replace(str(temp_merged) + ".csi", str(self.variants_bcf) + ".csi")

            # Log update details
            input_md5 = compute_md5(self.input_file)
            self.log_message(f"Added new file: {self.input_file}")
            self.log_message(f"Input file MD5: {input_md5}")
            self.log_message("Database statistics after update:")
            for key, value in stats.items():
                self.log_message(f"{key}: {value}")
            self.log_message(f"Processing completed in {duration.total_seconds():.2f} seconds")

        except subprocess.CalledProcessError as e:
            # Clean up temporary files on error
            temp_merged.unlink(missing_ok=True)
            Path(str(temp_merged) + ".csi").unlink(missing_ok=True)
            raise RuntimeError(f"Failed to add variants: {e}")