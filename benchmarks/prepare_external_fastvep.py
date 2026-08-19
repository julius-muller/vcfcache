#!/usr/bin/env python3
"""Build fastVEP caches from the exact VEP external-campaign blueprints.

The output is a strategy manifest accepted by ``run_external_cohort.py
prepare --tool fastvep``.  Cache membership is inherited from the frozen VEP
campaign, while every annotation is rebuilt with the supplied fastVEP recipe.
"""

from __future__ import annotations

import argparse
import json
import shlex
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from benchmarks.run_cohort import sha256sum
from benchmarks.run_pilot import write_json_atomic

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_VCFCACHE = REPO_ROOT / ".venv/bin/vcfcache"
DEFAULT_CACHE_NAME = "cache-fastvep-publication"


def _command_version(command: list[str]) -> str:
    completed = subprocess.run(command, capture_output=True, text=True, check=True)
    return (completed.stdout or completed.stderr).strip().splitlines()[0]


def _git_commit(path: Path) -> str | None:
    """Return the commit owning a tool path when it belongs to a Git checkout."""
    completed = subprocess.run(
        ["git", "-C", str(path.parent), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def _fingerprint(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Required immutable asset is missing: {path}")
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256sum(path),
    }


def freeze_toolchain(vcfcache: Path, params: Path) -> dict[str, Any]:
    """Freeze executables and every annotation asset referenced by params."""
    values = yaml.safe_load(params.read_text()) or {}
    if not isinstance(values, dict):
        raise RuntimeError(f"Params must be a mapping: {params}")
    fastvep_tokens = shlex.split(str(values.get("annotation_tool_cmd", "")))
    if not fastvep_tokens:
        raise RuntimeError(f"annotation_tool_cmd is absent from {params}")
    fastvep_resolved = shutil.which(fastvep_tokens[0]) or fastvep_tokens[0]
    fastvep = Path(fastvep_resolved).resolve()
    bcftools_tokens = shlex.split(str(values.get("bcftools_cmd", "bcftools")))
    bcftools_resolved = shutil.which(bcftools_tokens[0])
    if not bcftools_resolved:
        raise FileNotFoundError(f"bcftools is unavailable: {bcftools_tokens[0]}")
    bcftools = Path(bcftools_resolved).resolve()

    assets: dict[str, dict[str, Any]] = {}
    for name in ("fasta", "transcript_cache"):
        value = values.get(name)
        if not value:
            raise RuntimeError(f"Required {name} is absent from {params}")
        asset = Path(str(value)).resolve()
        assets[name] = _fingerprint(asset)
        if name == "fasta":
            for suffix in (".fai", ".gzi"):
                index = Path(f"{asset}{suffix}")
                if index.exists():
                    assets[f"fasta{suffix}"] = _fingerprint(index)

    def find_databases(value: Any, prefix: str = "params") -> None:
        if isinstance(value, dict):
            for key, item in value.items():
                find_databases(item, f"{prefix}.{key}")
        elif isinstance(value, list):
            for index, item in enumerate(value):
                find_databases(item, f"{prefix}[{index}]")
        elif isinstance(value, str):
            for token in shlex.split(value):
                if token.endswith((".osa", ".osa2")):
                    assets[prefix] = _fingerprint(Path(token).resolve())

    find_databases(values)
    vcfcache = vcfcache.resolve()
    return {
        "fastvep": {
            **_fingerprint(fastvep),
            "version": _command_version([str(fastvep), "--version"]),
            "git_commit": _git_commit(fastvep),
        },
        "vcfcache": {
            **_fingerprint(vcfcache),
            "version": _command_version([str(vcfcache), "--version"]),
            "git_commit": _git_commit(vcfcache),
        },
        "bcftools": {
            **_fingerprint(bcftools),
            "version": _command_version([str(bcftools), "--version"]),
        },
        "assets": assets,
        "params": _fingerprint(params),
    }


def source_blueprint(cache_dir: Path) -> Path:
    """Resolve the immutable blueprint belonging to a VCFcache cache."""
    path = cache_dir.parents[1] / "blueprint/vcfcache.bcf"
    if not path.exists() or not Path(f"{path}.csi").exists():
        raise FileNotFoundError(f"Source blueprint or index is missing: {path}")
    return path


def run(command: list[str | Path], log: Path) -> None:
    """Run one cache-construction command with an auditable log."""
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a") as handle:
        handle.write("$ " + " ".join(map(str, command)) + "\n")
        handle.flush()
        subprocess.run(
            list(map(str, command)),
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
            text=True,
        )


def published(path: Path, output_root: Path, published_root: Path) -> str:
    """Translate preparation paths to the worker-visible shared mount."""
    return str(published_root / path.relative_to(output_root))


def build_one(
    *,
    vcfcache: Path,
    output_root: Path,
    published_root: Path,
    key: str,
    source: dict[str, Any],
    recipe: Path,
    params: Path,
    cache_name: str,
    toolchain: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Rebuild one cache with unchanged membership and a fastVEP recipe."""
    assembly = source["assembly"]
    blueprint = source_blueprint(Path(source["cache_dir"]))
    database = output_root / "databases" / assembly / key
    cache = database / "cache" / cache_name
    output = cache / "vcfcache_annotated.bcf"
    provenance_path = database / "fastvep_cache_provenance.json"
    output_complete = output.exists() and Path(f"{output}.csi").exists()
    original: dict[str, Any] | None = None
    if output_complete:
        if not provenance_path.exists():
            raise RuntimeError(
                f"Refusing to invent build time for an existing cache: {cache}"
            )
        original = json.loads(provenance_path.read_text())
        if not original.get("complete") or not isinstance(
            original.get("wall_seconds"), (int, float)
        ):
            raise RuntimeError(f"Existing provenance is incomplete: {provenance_path}")
        if original.get("cache_sha256") != sha256sum(output):
            raise RuntimeError(f"Existing cache changed: {output}")
    started_at = datetime.now(timezone.utc)
    started = time.monotonic()
    if not (database / "blueprint/vcfcache.bcf").exists():
        run(
            [
                vcfcache,
                "blueprint-init",
                "--vcf",
                blueprint,
                "--output",
                database,
            ],
            output_root / f"logs/{key}.blueprint.log",
        )
    if not output_complete:
        run(
            [
                vcfcache,
                "cache-build",
                "--db",
                database,
                "--name",
                cache_name,
                "--anno-config",
                recipe,
                "--params",
                params,
            ],
            output_root / f"logs/{key}.cache.log",
        )
    if not output.exists() or not Path(f"{output}.csi").exists():
        raise RuntimeError(f"fastVEP cache build did not complete: {cache}")
    build_identity = {
        "tool": "fastvep",
        "kind": "fastvep_reannotation_of_frozen_blueprint",
        "assembly": assembly,
        "strategy": source["name"],
        "source_blueprint": str(blueprint),
        "source_blueprint_sha256": sha256sum(blueprint),
        "annotation_yaml_sha256": sha256sum(cache / "annotation.yaml"),
        "params_yaml_sha256": sha256sum(cache / "params.snapshot.yaml"),
        "cache_sha256": sha256sum(output),
        "complete": True,
    }
    if original is None:
        provenance = {
            **build_identity,
            "started_at": started_at.isoformat(),
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": time.monotonic() - started,
            "reused": False,
            "toolchain": toolchain,
        }
    else:
        for key, value in build_identity.items():
            if key in original and original[key] != value:
                raise RuntimeError(
                    f"Existing provenance mismatch at {key}: "
                    f"{original[key]!r} != {value!r}"
                )
        provenance = {
            **original,
            **build_identity,
            "reused": True,
            "last_verified_at": datetime.now(timezone.utc).isoformat(),
            "toolchain": toolchain,
        }
    write_json_atomic(provenance_path, provenance)
    expected = {
        key: provenance[key]
        for key in ("tool", "kind", "assembly", "strategy", "complete")
    }
    source_kind = source.get("kind", "unknown")
    local_kind = (
        "locally_built_fastvep_from_bundled_blueprint"
        if source_kind == "bundled_zenodo"
        else "locally_built_fastvep_from_cohort_blueprint"
    )
    strategy = {
        **{name: value for name, value in source.items() if name != "cache_dir"},
        "kind": local_kind,
        "cache_dir": published(cache, output_root, published_root),
        "controller_cache_dir": str(cache),
        "annotation_yaml_sha256": provenance["annotation_yaml_sha256"],
        "provenance_path": published(provenance_path, output_root, published_root),
        "controller_provenance_path": str(provenance_path),
        "provenance_expected": expected,
        "source_cache_dir": source["cache_dir"],
        "source_strategy_kind": source_kind,
        "source_alias": source.get("alias", ""),
        "source_doi": source.get("doi", ""),
        "alias": "",
        "doi": "",
        "source_blueprint_sha256": provenance["source_blueprint_sha256"],
        "build_wall_seconds": provenance["wall_seconds"],
    }
    return strategy, provenance


def prepare(args: argparse.Namespace) -> Path:
    """Materialize all seven tool-specific caches and freeze their manifest."""
    source = json.loads(args.vep_strategies.read_text())
    cohort_assemblies = source["cohort_assemblies"]
    if set(cohort_assemblies.values()) != {"GRCh37", "GRCh38"}:
        raise RuntimeError("Exact VEP parity requires both GRCh37 and GRCh38")
    configs = {
        "GRCh37": (args.recipe_grch37, args.params_grch37),
        "GRCh38": (args.recipe_grch38, args.params_grch38),
    }
    for assembly, paths in configs.items():
        for path in paths:
            if not path.exists():
                raise FileNotFoundError(f"Missing {assembly} fastVEP config: {path}")
    if not args.vcfcache.exists():
        raise FileNotFoundError(f"VCFcache executable is missing: {args.vcfcache}")
    args.output_root.mkdir(parents=True, exist_ok=True)

    built: dict[tuple[str, str], dict[str, Any]] = {}
    provenance: dict[str, Any] = {}
    for assembly, strategies in source["bundled_strategies_by_assembly"].items():
        recipe, params = configs[assembly]
        toolchain = freeze_toolchain(args.vcfcache, params)
        for item in strategies:
            key = f"bundled-{assembly}-{item['name']}"
            built[(assembly, item["name"])], provenance[key] = build_one(
                vcfcache=args.vcfcache,
                output_root=args.output_root,
                published_root=args.published_root,
                key=key,
                source=item,
                recipe=recipe,
                params=params,
                cache_name=args.cache_name,
                toolchain=toolchain,
            )

    custom: dict[str, dict[str, Any]] = {}
    for cohort, item in source["cohort_strategies"].items():
        assembly = cohort_assemblies[cohort]
        recipe, params = configs[assembly]
        toolchain = freeze_toolchain(args.vcfcache, params)
        key = f"cohort-{cohort}"
        custom[cohort], provenance[key] = build_one(
            vcfcache=args.vcfcache,
            output_root=args.output_root,
            published_root=args.published_root,
            key=key,
            source=item,
            recipe=recipe,
            params=params,
            cache_name=args.cache_name,
            toolchain=toolchain,
        )

    runtime: dict[str, dict[str, str]] = {}
    for assembly, (_recipe, params) in configs.items():
        destination = args.output_root / f"config/runtime_params_{assembly}.yaml"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(params, destination)
        runtime[assembly] = {
            "path": published(destination, args.output_root, args.published_root),
            "controller_path": str(destination),
            "sha256": sha256sum(destination),
            "source": str(params),
            "source_sha256": sha256sum(params),
        }

    document = {
        "tool": "fastvep",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_vep_strategies": str(args.vep_strategies),
        "source_vep_strategies_sha256": sha256sum(args.vep_strategies),
        "membership_identical_to_vep": True,
        "cohort_assemblies": cohort_assemblies,
        "bundled_strategies_by_assembly": {
            assembly: [built[(assembly, item["name"])] for item in strategies]
            for assembly, strategies in source["bundled_strategies_by_assembly"].items()
        },
        "cohort_strategies": custom,
        "runtime_params_by_assembly": runtime,
        "cache_builds": provenance,
    }
    output = args.output_root / "fastvep_strategies.json"
    write_json_atomic(output, document)
    print(output)
    return output


def parser() -> argparse.ArgumentParser:
    """Build the fastVEP cache-preparation CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--vep-strategies", type=Path, required=True)
    result.add_argument("--output-root", type=Path, required=True)
    result.add_argument(
        "--published-root",
        type=Path,
        help="Worker-visible equivalent of output-root (defaults to output-root)",
    )
    result.add_argument("--recipe-grch37", type=Path, required=True)
    result.add_argument("--params-grch37", type=Path, required=True)
    result.add_argument("--recipe-grch38", type=Path, required=True)
    result.add_argument("--params-grch38", type=Path, required=True)
    result.add_argument("--vcfcache", type=Path, default=DEFAULT_VCFCACHE)
    result.add_argument("--cache-name", default=DEFAULT_CACHE_NAME)
    return result


def main() -> None:
    """Prepare the annotator-specific cache bundle."""
    args = parser().parse_args()
    for name, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, name, value.expanduser().resolve())
    if args.published_root is None:
        args.published_root = args.output_root
    prepare(args)


if __name__ == "__main__":
    main()
