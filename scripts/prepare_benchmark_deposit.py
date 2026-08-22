#!/usr/bin/env python3
"""Assemble the benchmark archive for Zenodo and report what it will contain.

The manuscript cites a benchmark deposit that must let a reader re-derive every
reported number. This collects the frozen source tables, manifests, analysis
code and figure sources into one staging tree, writes a checksum manifest, and
prints the totals so the deposit can be reviewed before anything is uploaded.

Reserve the DOI on a Zenodo draft first, paste it into the manuscript, then
publish: the DOI has to exist before the text that cites it is final.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
from pathlib import Path

# Paths are relative to the repository root. Missing entries are reported rather
# than skipped silently, because a partial deposit is worse than a failed one.
CONTENTS = [
    ("figure source tables", "drafts/gigascience_technical_note/figures/source_data"),
    ("benchmark source tables", "benchmarks/figures/source_data"),
    ("claim ledger", "drafts/gigascience_technical_note/claims.tsv"),
    ("figure provenance", "drafts/gigascience_technical_note/figures/PROVENANCE.tsv"),
    ("assay campaign metrics", "benchmarks/external_assay_v3/source_data"),
    ("fastVEP pilot metrics", "benchmarks/fastvep_pilot/source_data"),
    ("analysis code", "benchmarks/analyze_breakeven.py"),
    ("analysis code", "benchmarks/analyze_richness_ladder.py"),
    ("analysis code", "benchmarks/analyze_assay_benchmark.py"),
    ("analysis code", "benchmarks/analyze_external_benchmark.py"),
    ("materials and methods", "benchmarks/MATERIALS_AND_METHODS.md"),
    ("storage layout", "benchmarks/docs/STORAGE_LAYOUT.md"),
    ("source provenance", "benchmarks/SOURCE_PROVENANCE.md"),
    ("runnable example data", "vcfcache/demo_data"),
]


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of one file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stage(root: Path, dest: Path) -> tuple[list[Path], list[str]]:
    """Copy every listed item into the staging tree, reporting what is absent."""
    copied: list[Path] = []
    missing: list[str] = []
    for label, rel in CONTENTS:
        src = root / rel
        if not src.exists():
            missing.append(f"{label}: {rel}")
            continue
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if src.is_dir():
            shutil.copytree(src, target, dirs_exist_ok=True)
            copied.extend(p for p in target.rglob("*") if p.is_file())
        else:
            shutil.copy2(src, target)
            copied.append(target)
    return copied, missing


def main() -> None:
    """Stage the deposit and write its checksum manifest."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--dest", type=Path, required=True)
    args = parser.parse_args()

    if args.dest.exists():
        shutil.rmtree(args.dest)
    args.dest.mkdir(parents=True)

    copied, missing = stage(args.root, args.dest)
    manifest = args.dest / "CHECKSUMS.sha256"
    total = 0
    with manifest.open("w") as handle:
        for path in sorted(copied):
            total += path.stat().st_size
            handle.write(f"{sha256(path)}  {path.relative_to(args.dest)}\n")

    print(f"staged {len(copied)} files, {total / 1e6:.1f} MB -> {args.dest}")
    print(f"manifest: {manifest}")
    if missing:
        print(f"\nMISSING ({len(missing)}) — resolve before depositing:")
        for item in missing:
            print(f"  {item}")
    else:
        print("\nall listed contents present")


if __name__ == "__main__":
    main()
