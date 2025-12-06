import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def run_cmd(cmd, cwd=None):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    assert result.returncode == 0, f"Command failed: {cmd}\nSTDOUT:{result.stdout}\nSTDERR:{result.stderr}"
    return result.stdout


def test_full_integration_annotation():
    """
    Integration test: download cache via alias + manifest, run annotate, and validate output exists.

    This test is enabled by default and relies on network access to Zenodo.
    It uses the public cache alias GRCh38-af010-vep115.2_basic.
    """

    alias = "GRCh38-af010-vep115.2_basic"
    manifest = Path(__file__).resolve().parent.parent / "public_caches.yaml"
    sample_vcf = Path(__file__).resolve().parent / "data" / "nodata" / "sample4.bcf"

    outdir = Path(tempfile.mkdtemp(prefix="vcfstash_integration_"))

    try:
        cmd = (
            f"vcfstash annotate -a {alias} "
            f"--vcf {sample_vcf} "
            f"--output {outdir} "
            f"--manifest {manifest} "
            f"--force "
        )
        run_cmd(cmd)

        # Validate output exists
        produced = outdir / f"{sample_vcf.stem}_vst.bcf"
        assert produced.exists(), f"Annotated BCF missing: {produced}"

        # bcftools sanity check
        run_cmd(f"bcftools view -h {produced}")

    finally:
        shutil.rmtree(outdir, ignore_errors=True)

