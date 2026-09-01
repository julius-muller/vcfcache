#!/usr/bin/env python3
"""Freeze checksums and renderer provenance for the selected figure set."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

DRAFT_DIR = Path(__file__).resolve().parents[1]
PROVENANCE = DRAFT_DIR / "figures/PROVENANCE.tsv"
FORMATS = ("svg", "pdf", "png")


STORY = "benchmarks/figures/output/final/2026-08-27T225600Z-story-campaign-final-v4"

FIGURES = (
    ("figures/main_v4/figure1_workflow_and_value_v4", "main",
     f"{STORY}/manuscript/figure_workflow_and_value", True,
     "recipe-specific lookup workflow and value drivers"),
    ("figures/main_v4/figure2_empirical_time_and_throughput_v4", "main",
     f"{STORY}/manuscript/figure_empirical_time_and_throughput", True,
     "matched native-annotator runtime remaining and serial throughput"),
    ("figures/main_v4/figure3_runtime_rule_v4", "main",
     f"{STORY}/manuscript/figure_runtime_rule", True,
     "empirical cache-coverage rule and controlled VEP hours returned"),
    ("figures/main_v4/figure4_cohort_economics_v4", "main",
     f"{STORY}/manuscript/figure_cohort_economics", True,
     "three-genome cohort-cache amortization and compute days returned"),
    ("figures/main_v4/panels/figure1_A_workflow", "main_panel",
     f"{STORY}/manuscript/panels/figure_workflow_A", True,
     "standalone Figure 1A"),
    ("figures/main_v4/panels/figure1_B_value_drivers", "main_panel",
     f"{STORY}/manuscript/panels/figure_workflow_B", True,
     "standalone Figure 1B"),
    ("figures/main_v4/panels/figure2_A_runtime_remaining", "main_panel",
     f"{STORY}/manuscript/panels/figure_empirical_A_remaining", True,
     "standalone Figure 2A"),
    ("figures/main_v4/panels/figure2_B_throughput", "main_panel",
     f"{STORY}/manuscript/panels/figure_empirical_B_throughput", True,
     "standalone Figure 2B"),
    ("figures/main_v4/panels/figure3_A_empirical_rule", "main_panel",
     f"{STORY}/manuscript/panels/figure_runtime_rule_A_empirical", True,
     "standalone Figure 3A"),
    ("figures/main_v4/panels/figure3_B_pipeline_cost", "main_panel",
     f"{STORY}/manuscript/panels/figure_runtime_rule_B_complexity", True,
     "standalone Figure 3B"),
    ("figures/main_v4/panels/figure4_A_amortization", "main_panel",
     f"{STORY}/manuscript/panels/figure_economics_A_amortization", True,
     "standalone Figure 4A"),
    ("figures/main_v4/panels/figure4_B_days_returned", "main_panel",
     f"{STORY}/manuscript/panels/figure_economics_B_days_returned", True,
     "standalone Figure 4B"),
    ("supplement/figures/supplementary_figure1_external_wgs_v2", "supplement",
     "drafts/gigascience_technical_note/scripts/render_publication_figures.R", False,
     "cache coverage and wall time returned in 52 independent WGS samples"),
    ("supplement/figures/panels/supplementary_figure1_A_cache_coverage",
     "supplement_panel",
     "drafts/gigascience_technical_note/scripts/render_publication_figures.R", False,
     "standalone Supplementary Figure 1A"),
    ("supplement/figures/panels/supplementary_figure1_B_time_saved", "supplement_panel",
     "drafts/gigascience_technical_note/scripts/render_publication_figures.R", False,
     "standalone Supplementary Figure 1B"),
    ("supplement/figures/supplementary_figure2_fastvep_cpu_complexity", "supplement",
     "drafts/gigascience_technical_note/scripts/render_fastvep_cpu_complexity.R", False,
     "process-wide CPU fastVEP workload-complexity interaction"),
    ("supplement/figures/panels/supplementary_figure2_A_runtime", "supplement_panel",
     "drafts/gigascience_technical_note/scripts/render_fastvep_cpu_complexity.R", False,
     "standalone Supplementary Figure 2A"),
    ("supplement/figures/panels/supplementary_figure2_B_speedup", "supplement_panel",
     "drafts/gigascience_technical_note/scripts/render_fastvep_cpu_complexity.R", False,
     "standalone Supplementary Figure 2B"),
    ("supplement/figures/v4/supplementary_figure3_annotation_load_fidelity_v4",
     "supplement",
     f"{STORY}/supplement/figure_annotation_load_and_fidelity", True,
     "paired SnpEff annotation load and final semantic audit"),
    ("supplement/figures/v4/panels/supplementary_figure3_A_snpeff_load",
     "supplement_panel",
     f"{STORY}/supplement/panels/figure_load_A_snpeff", True,
     "standalone Supplementary Figure 3A"),
    ("supplement/figures/v4/panels/supplementary_figure3_B_fidelity",
     "supplement_panel",
     f"{STORY}/supplement/panels/figure_fidelity_B_audit", True,
     "standalone Supplementary Figure 3B"),
    ("figures/candidates_v4/headliner_one_day_of_compute", "candidate",
     f"{STORY}/repository/headliner_one_day_of_compute", True,
     "repository serial-throughput headline; not cited"),
    ("figures/candidates_v4/headliner_wgs_clock", "candidate",
     f"{STORY}/repository/headliner_wgs_clock", True,
     "repository WGS direct-to-cached clock; not cited"),
    ("figures/candidates_v4/cohort_100_time_returned", "candidate",
     f"{STORY}/alternatives/cohort_100_time_returned", True,
     "presentation-style 100-genome cohort value graphic; not cited"),
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
    for prefix, role, source, source_is_asset, note in FIGURES:
        for suffix in FORMATS:
            relative = Path(f"{prefix}.{suffix}")
            path = DRAFT_DIR / relative
            if not path.is_file():
                raise SystemExit(f"Selected figure asset is missing: {path}")
            rows.append(
                {
                    "target": relative.as_posix(),
                    "role": role,
                    "source": f"{source}.{suffix}" if source_is_asset else source,
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
