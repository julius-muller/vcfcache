from pathlib import Path

from benchmarks.fastvep_pilot.run import (
    Config,
    canonical_info,
    vcfcache_command,
)
from benchmarks.prepare_fastvep_node import write_configs, write_fasta_index


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


def test_fastvep_fasta_index_does_not_require_samtools(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_bytes(b">chr1 description\nACGT\nTG\n>chr2\nAAA\n")
    index = tmp_path / "reference.fa.fai"
    write_fasta_index(fasta, index)
    assert index.read_text().splitlines() == [
        "chr1\t6\t18\t4\t5",
        "chr2\t3\t32\t3\t4",
    ]


def test_fastvep_publication_recipe_omits_known_noop_flags(tmp_path: Path) -> None:
    recipe, _params = write_configs(
        tmp_path,
        "GRCh38",
        Path("/tools/fastvep"),
        Path("/data/transcripts.cache"),
        Path("/data/reference.fa"),
        8,
    )
    text = recipe.read_text()
    assert "--hgvs --no-progress" in text
    assert "--everything" not in text
    assert "--symbol" not in text
    assert "--canonical" not in text
