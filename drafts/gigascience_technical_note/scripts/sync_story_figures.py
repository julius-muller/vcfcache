#!/usr/bin/env python3
"""Copy the frozen story-driven figure snapshot into the manuscript workspace."""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

DRAFT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = DRAFT_DIR.parents[1]
SNAPSHOT = (
    REPO_ROOT
    / "benchmarks/figures/output/final/"
    "2026-08-27T225600Z-story-campaign-final-v4"
)
FORMATS = ("svg", "pdf", "png")

PREFIXES = {
    "manuscript/figure_workflow_and_value":
        "figures/main_v4/figure1_workflow_and_value_v4",
    "manuscript/figure_empirical_time_and_throughput":
        "figures/main_v4/figure2_empirical_time_and_throughput_v4",
    "manuscript/figure_runtime_rule":
        "figures/main_v4/figure3_runtime_rule_v4",
    "manuscript/figure_cohort_economics":
        "figures/main_v4/figure4_cohort_economics_v4",
    "manuscript/panels/figure_workflow_A":
        "figures/main_v4/panels/figure1_A_workflow",
    "manuscript/panels/figure_workflow_B":
        "figures/main_v4/panels/figure1_B_value_drivers",
    "manuscript/panels/figure_empirical_A_remaining":
        "figures/main_v4/panels/figure2_A_runtime_remaining",
    "manuscript/panels/figure_empirical_B_throughput":
        "figures/main_v4/panels/figure2_B_throughput",
    "manuscript/panels/figure_runtime_rule_A_empirical":
        "figures/main_v4/panels/figure3_A_empirical_rule",
    "manuscript/panels/figure_runtime_rule_B_complexity":
        "figures/main_v4/panels/figure3_B_pipeline_cost",
    "manuscript/panels/figure_economics_A_amortization":
        "figures/main_v4/panels/figure4_A_amortization",
    "manuscript/panels/figure_economics_B_days_returned":
        "figures/main_v4/panels/figure4_B_days_returned",
    "supplement/figure_annotation_load_and_fidelity":
        "supplement/figures/v4/supplementary_figure3_annotation_load_fidelity_v4",
    "supplement/panels/figure_load_A_snpeff":
        "supplement/figures/v4/panels/supplementary_figure3_A_snpeff_load",
    "supplement/panels/figure_fidelity_B_audit":
        "supplement/figures/v4/panels/supplementary_figure3_B_fidelity",
    "repository/headliner_one_day_of_compute":
        "figures/candidates_v4/headliner_one_day_of_compute",
    "repository/headliner_wgs_clock":
        "figures/candidates_v4/headliner_wgs_clock",
    "alternatives/cohort_100_time_returned":
        "figures/candidates_v4/cohort_100_time_returned",
}


def sha256(path: Path) -> str:
    """Return the SHA-256 checksum of one file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_checked(source: Path, target: Path) -> None:
    """Copy one frozen asset and verify that its content is unchanged."""
    if not source.is_file():
        raise SystemExit(f"Frozen source is missing: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    if sha256(source) != sha256(target):
        raise SystemExit(f"Copied figure checksum differs: {target}")


def main() -> None:
    """Synchronize composites, panels, source tables, and snapshot metadata."""
    if not SNAPSHOT.is_dir():
        raise SystemExit(f"Frozen story snapshot is missing: {SNAPSHOT}")

    for source_prefix, target_prefix in PREFIXES.items():
        for suffix in FORMATS:
            copy_checked(
                SNAPSHOT / f"{source_prefix}.{suffix}",
                DRAFT_DIR / f"{target_prefix}.{suffix}",
            )

    source_target = DRAFT_DIR / "figures/source_data/story_v4"
    source_target.mkdir(parents=True, exist_ok=True)
    for source in sorted((SNAPSHOT / "source").glob("*.tsv")):
        copy_checked(source, source_target / source.name)

    for name in (
        "FIGURE_INDEX.md",
        "OUTPUT_MANIFEST.tsv",
        "SOURCE_PROVENANCE.tsv",
        "R_SESSION_INFO.txt",
    ):
        copy_checked(SNAPSHOT / name, DRAFT_DIR / "figures/main_v4" / name)

    print(f"Synchronized {len(PREFIXES)} figure prefixes from {SNAPSHOT.name}")


if __name__ == "__main__":
    main()
