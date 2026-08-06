#!/usr/bin/env python3
"""Build fastVEP caches from the exact VEP external-campaign blueprints.

The output is a strategy manifest accepted by ``run_external_cohort.py
prepare --tool fastvep``.  Cache membership is inherited from the frozen VEP
campaign, while every annotation is rebuilt with the supplied fastVEP recipe.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from benchmarks.run_cohort import sha256sum
from benchmarks.run_pilot import write_json_atomic

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_VCFCACHE = REPO_ROOT / ".venv/bin/vcfcache"
DEFAULT_CACHE_NAME = "cache-fastvep-publication"


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
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Rebuild one cache with unchanged membership and a fastVEP recipe."""
    assembly = source["assembly"]
    blueprint = source_blueprint(Path(source["cache_dir"]))
    database = output_root / "databases" / assembly / key
    cache = database / "cache" / cache_name
    output = cache / "vcfcache_annotated.bcf"
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
    if not output.exists() or not Path(f"{output}.csi").exists():
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
    provenance = {
        "tool": "fastvep",
        "kind": "fastvep_reannotation_of_frozen_blueprint",
        "assembly": assembly,
        "strategy": source["name"],
        "source_blueprint": str(blueprint),
        "source_blueprint_sha256": sha256sum(blueprint),
        "annotation_yaml_sha256": sha256sum(cache / "annotation.yaml"),
        "params_yaml_sha256": sha256sum(cache / "params.snapshot.yaml"),
        "cache_sha256": sha256sum(output),
        "started_at": started_at.isoformat(),
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "wall_seconds": time.monotonic() - started,
        "complete": True,
    }
    provenance_path = database / "fastvep_cache_provenance.json"
    write_json_atomic(provenance_path, provenance)
    expected = {
        key: provenance[key]
        for key in ("tool", "kind", "assembly", "strategy", "complete")
    }
    strategy = {
        **{name: value for name, value in source.items() if name != "cache_dir"},
        "cache_dir": published(cache, output_root, published_root),
        "controller_cache_dir": str(cache),
        "annotation_yaml_sha256": provenance["annotation_yaml_sha256"],
        "provenance_path": published(provenance_path, output_root, published_root),
        "controller_provenance_path": str(provenance_path),
        "provenance_expected": expected,
        "source_cache_dir": source["cache_dir"],
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
            )

    custom: dict[str, dict[str, Any]] = {}
    for cohort, item in source["cohort_strategies"].items():
        assembly = cohort_assemblies[cohort]
        recipe, params = configs[assembly]
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
