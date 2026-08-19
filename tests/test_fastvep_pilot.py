import json
from pathlib import Path

from benchmarks.fastvep_pilot.core_scaling import (
    copy_required_cache_structure,
    write_controlled_cache_metadata,
)
from benchmarks.fastvep_pilot.run import (
    Config,
    canonical_info,
    vcfcache_command,
)
from benchmarks.prepare_external_fastvep import build_one
from benchmarks.prepare_fastvep_node import write_configs, write_fasta_index
from benchmarks.run_cohort import sha256sum


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


def test_controlled_fastvep_cache_writes_validation_metadata(tmp_path: Path) -> None:
    cache = tmp_path / "cache"
    cache.mkdir()
    input_vcf = tmp_path / "independent.vcf.gz"
    write_controlled_cache_metadata(cache, input_vcf, 4_316_659)
    metadata = json.loads((cache / "cache.json").read_text())
    assert metadata == {
        "construction": "POS modulo 100 < 90",
        "input": str(input_vcf),
        "records": 4_316_659,
        "target_hit_rate": 0.9,
    }


def test_controlled_fastvep_cache_copies_frozen_root_contract(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source"
    (source / "blueprint").mkdir(parents=True)
    (source / "blueprint/vcfcache.bcf").write_bytes(b"blueprint")
    (source / "workflow").mkdir()
    (source / "workflow/init.yaml").write_text("genome_build: GRCh38\n")
    destination = tmp_path / "destination"
    copy_required_cache_structure(source, destination)
    assert (destination / "blueprint/vcfcache.bcf").read_bytes() == b"blueprint"
    assert (destination / "workflow/init.yaml").read_text() == (
        "genome_build: GRCh38\n"
    )


def test_existing_fastvep_cache_preserves_build_time_and_local_provenance(tmp_path):
    source_cache = tmp_path / "source/database/cache/vep"
    source_cache.mkdir(parents=True)
    source_blueprint = tmp_path / "source/database/blueprint/vcfcache.bcf"
    source_blueprint.parent.mkdir(parents=True)
    source_blueprint.write_bytes(b"membership")
    Path(f"{source_blueprint}.csi").write_bytes(b"index")

    output_root = tmp_path / "output"
    database = output_root / "databases/GRCh38/bundled-GRCh38-gnomad_af_0.01"
    local_blueprint = database / "blueprint/vcfcache.bcf"
    local_blueprint.parent.mkdir(parents=True)
    local_blueprint.write_bytes(b"membership")
    cache = database / "cache/cache-fastvep"
    cache.mkdir(parents=True)
    output = cache / "vcfcache_annotated.bcf"
    output.write_bytes(b"annotated")
    Path(f"{output}.csi").write_bytes(b"index")
    (cache / "annotation.yaml").write_text("genome_build: GRCh38\n")
    (cache / "params.snapshot.yaml").write_text("threads: 8\n")
    provenance_path = database / "fastvep_cache_provenance.json"
    provenance_path.write_text(
        __import__("json").dumps(
            {
                "tool": "fastvep",
                "kind": "fastvep_reannotation_of_frozen_blueprint",
                "assembly": "GRCh38",
                "strategy": "gnomad_af_0.01",
                "source_blueprint": str(source_blueprint),
                "source_blueprint_sha256": sha256sum(source_blueprint),
                "annotation_yaml_sha256": sha256sum(cache / "annotation.yaml"),
                "params_yaml_sha256": sha256sum(cache / "params.snapshot.yaml"),
                "cache_sha256": sha256sum(output),
                "started_at": "2026-08-01T10:00:00+00:00",
                "completed_at": "2026-08-01T10:10:00+00:00",
                "wall_seconds": 600.0,
                "complete": True,
            }
        )
    )
    recipe = tmp_path / "annotation.yaml"
    params = tmp_path / "params.yaml"
    recipe.write_text("recipe\n")
    params.write_text("params\n")
    strategy, provenance = build_one(
        vcfcache=tmp_path / "vcfcache",
        output_root=output_root,
        published_root=output_root,
        key="bundled-GRCh38-gnomad_af_0.01",
        source={
            "name": "gnomad_af_0.01",
            "assembly": "GRCh38",
            "kind": "bundled_zenodo",
            "alias": "cache-original",
            "doi": "10.example/original",
            "cache_dir": str(source_cache),
        },
        recipe=recipe,
        params=params,
        cache_name="cache-fastvep",
        toolchain={"frozen": True},
    )
    assert provenance["wall_seconds"] == 600.0
    assert provenance["started_at"] == "2026-08-01T10:00:00+00:00"
    assert provenance["reused"] is True
    assert strategy["build_wall_seconds"] == 600.0
    assert strategy["kind"] == "locally_built_fastvep_from_bundled_blueprint"
    assert strategy["source_strategy_kind"] == "bundled_zenodo"
    assert strategy["source_alias"] == "cache-original"
    assert strategy["alias"] == ""
