#!/usr/bin/env python3
"""Freeze checksums and renderer provenance for the selected figure set."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

DRAFT_DIR = Path(__file__).resolve().parents[1]
PROVENANCE = DRAFT_DIR / "figures/PROVENANCE.tsv"
FORMATS = ("svg", "pdf", "png")


FIGURES = (
    (
        "figures/main/figure1_workflow_v3",
        "main",
        "drafts/gigascience_technical_note/scripts/render_publication_figures_v3.R",
        "Zenodo, cohort, and per-sample recipe-specific cache workflow",
    ),
    (
        "figures/main/figure2_assay_annotator_v3",
        "main",
        "drafts/gigascience_technical_note/scripts/render_publication_figures_v3.R",
        "matched Panel, WES, and WGS use-case boundary for VEP and fastVEP",
    ),
    (
        "figures/main/figure3_pipeline_complexity_v3",
        "main",
        "drafts/gigascience_technical_note/scripts/render_publication_figures_v3.R",
        "absolute time returned and hit-rate speedup ceilings",
    ),
    (
        "supplement/figures/supplementary_figure1_external_wgs_v2",
        "supplement",
        "drafts/gigascience_technical_note/scripts/render_publication_figures.R",
        "cache coverage and wall time returned in 52 independent WGS samples",
    ),
    (
        "supplement/figures/supplementary_figure2_fastvep_cpu_complexity",
        "supplement",
        "drafts/gigascience_technical_note/scripts/render_fastvep_cpu_complexity.R",
        "enforced-CPU fastVEP workload-complexity interaction in held-out KPGP WGS",
    ),
    (
        "supplement/figures/supplementary_figure3_matched_assays_v3",
        "supplement",
        "drafts/gigascience_technical_note/scripts/render_publication_figures_v3.R",
        "sample-level cache coverage and speedup in the matched assay extension",
    ),
    (
        "figures/candidates/headliner_before_after",
        "candidate",
        "benchmarks/figures/output/final/2026-08-09T223718Z-light-matched-final/repository/headliner_before_after",
        "simplified measured VEP and fastVEP headline; not cited",
    ),
    (
        "figures/candidates/headliner_workflow_doodle",
        "candidate",
        "benchmarks/figures/output/final/2026-08-09T223718Z-light-matched-final/repository/headliner_workflow_doodle",
        "repository-style workflow headline; not cited",
    ),
    (
        "figures/candidates/headliner_cache_coverage",
        "candidate",
        "benchmarks/figures/output/final/2026-08-12T223200Z-pipeline-complexity-dual-cache-final-v2/repository/headliner_cache_coverage",
        "simplified cache-depth headline; not cited",
    ),
)


def sha256(path: Path) -> str:
    """Calculate one figure asset digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    """Write the complete selected-figure provenance table."""
    rows: list[dict[str, str | int]] = []
    for prefix, role, source, note in FIGURES:
        for suffix in FORMATS:
            relative = Path(f"{prefix}.{suffix}")
            path = DRAFT_DIR / relative
            if not path.is_file():
                raise SystemExit(f"Selected figure asset is missing: {path}")
            rows.append(
                {
                    "target": relative.as_posix(),
                    "role": role,
                    "source": f"{source}.{suffix}" if role == "candidate" else source,
                    "sha256": sha256(path),
                    "bytes": path.stat().st_size,
                    "note": note,
                }
            )

    with PROVENANCE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("target", "role", "source", "sha256", "bytes", "note"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
