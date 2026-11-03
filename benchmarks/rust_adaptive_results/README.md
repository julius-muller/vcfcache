# Rust-adaptive manuscript figures (historical measured campaign)

This is a matching figure set for the **Rust-adaptive** implementation, not the
native/resident engine and not Python-adaptive timings with different labels.
No benchmarks were rerun. The original evidence comes from the frozen
2026-10-29 consolidation on `feature/performance-sprint-20261027`.

- **Main figures, supplements and repository headliner:** `final-20261103/figures-v1/`.
- **Prominent core versus annotation-rich vep-rs:**
  `interpretation-20261103/figures-v1/provisional/Figure4_veprs_core_and_annotation_rich`.
- **General sensitivity model:** `interpretation-20261103/figures-v1/model/Figure5_runtime_model_sensitivity`.
  It is the same hypothetical model used for the Python set, not fitted to Rust.
- **Data and RStudio tibbles:** `final-20261103/tables/` and
  `interpretation-20261103/tables/`.

The shared R/ggplot2 renderer uses identical filters, layouts, cohort-stratified
bootstrap procedure, axis semantics and output formats as the Python set. Each
multi-panel figure includes separate panel files. Assembly/recipe strata remain
separate; no 1000 Genomes cohort is plotted. Existing figures were not overwritten.

Median GRCh38 assay speedups, 12 samples per cell:

| Recipe | Panel | WES | WGS |
|---|---:|---:|---:|
| VEP everything | 2.59× | 2.46× | 10.23× |
| fastVEP clinical | 1.00× | 0.81× | 2.88× |
| SnpEff clinical | 1.27× | 2.65× | 10.81× |
| vep-rs core | 1.33× | 0.46× | 1.32× |
| vep-rs clinical | 1.21× | **1.93× timing only** | **3.18× timing only** |

The same 572 accepted, 24 clinical vep-rs timing-only and four excluded SnpEff
conditions are retained. Clinical CSQ mismatches are not excused to strengthen
the headline. The annotation-rich recipe genuinely has higher annotation cost,
but those 24 conditions cannot yet support validated speedup claims.

`adapt.py` checks the SHA256-bound original evidence archive, all 600 measured
elapsed times, CPU allocations and Rust helper identities. It also verifies that
all canonical conditions match the historical Rust controls bound into the Python
replay, including identical direct denominators, sample and cache hashes. It only
maps labels and memory units for the shared renderer; no timing is recomputed or
selected for favorable performance. Original exports and the evidence archive
remain in `final-20261103/source/original/`.

Comparisons between engines are historical, not a newly randomized A/B experiment;
date, node and filesystem-cache effects remain possible. Do not interpret the two
sets as evidence that every workload benefits from compiled code. For example,
core WES remains slower in both implementations.

The detailed captions and model specification are in
`../python_adaptive_results/interpretation-20261103/README.md` (Figures 4–5), and
`../python_adaptive_results/final-20261103/CAPTIONS.md` (Figures 1–3 and S1–S3).
Apply those descriptions to the Rust dataset here, not Python values. Figure 1
is the same unchanged conceptual workflow. The original Rust campaign's direct
controls were collected historically; producing this figure set launched zero jobs.

Reproduce:

```bash
python -m benchmarks.rust_adaptive_results.adapt
Rscript --vanilla benchmarks/python_adaptive_results/render.R \
  benchmarks/rust_adaptive_results/final-20261103 figures-v2
Rscript --vanilla benchmarks/python_adaptive_results/extend_analysis.R \
  benchmarks/rust_adaptive_results/final-20261103 \
  benchmarks/rust_adaptive_results/interpretation-20261103 figures-v2
```

Use new figure-directory names on every render. No default-engine change or
additional Slurm workload was introduced by this analysis.
