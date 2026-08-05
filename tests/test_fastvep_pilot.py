from pathlib import Path

from benchmarks.fastvep_pilot.run import (
    Config,
    canonical_info,
    vcfcache_command,
)


def test_canonical_info_sorts_tags_and_csq_entries() -> None:
    assert canonical_info("ZZ=2;CSQ=b,a;AA=1") == "AA=1;CSQ=a,b;ZZ=2"


def test_matched_commands_only_add_uncached_flag(tmp_path: Path) -> None:
    config = Config(
        root=tmp_path / "runtime",
        repo=tmp_path / "repo",
        input_vcf=tmp_path / "input.vcf.gz",
        fasta_gz=tmp_path / "reference.fa.gz",
        vep_sif=tmp_path / "vep.sif",
        threads=16,
    )
    cache = tmp_path / "cache"
    run_dir = tmp_path / "run"
    direct = [
        str(value) for value in vcfcache_command(config, cache, run_dir, uncached=True)
    ]
    cached = [
        str(value) for value in vcfcache_command(config, cache, run_dir, uncached=False)
    ]
    assert direct[:-1] == cached
    assert direct[-1] == "--uncached"
    assert "--skip-split-multiallelic" in direct
