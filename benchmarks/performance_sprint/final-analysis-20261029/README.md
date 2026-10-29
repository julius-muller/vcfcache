# Complete campaign analysis — 29 October 2026

Start with [the decision brief](DECISIONS.md), then the comprehensive
[analysis](ANALYSIS.md). All approved annotation jobs and collectors are finished.
No new annotation, cache changes, release defaults or manuscript edits were
performed as part of this consolidation.

## Key files

- [ANALYSIS.md](ANALYSIS.md): findings, uncertainty, limitations, bottlenecks,
  native/adaptive trade-offs and recommendations.
- [DECISIONS.md](DECISIONS.md): next decisions and explicit stop/go gates.
- [FIGURE_LEGENDS.md](FIGURE_LEGENDS.md): numbered captions and abbreviations.
- [source/SNAPSHOT.json](source/SNAPSHOT.json): source campaigns, hashes and
  terminal counts; [compressed evidence](source/terminal-evidence.tar.gz) retains
  selected original JSON/YAML files, including superseded attempts.
- [tables/cells.tsv](tables/cells.tsv): **600** canonical conditions.
- [tables/plot_views.tsv](tables/plot_views.tsv): **612** views, including 12
  aliases; use `measurement_id` to avoid double counting.
- [tables/manuscript_data.rds](tables/manuscript_data.rds): named, analysis-ready
  tibbles for RStudio, including raw rows, accepted/excluded assays, WGS
  strategies, phase timings, resources, provenance and optimization scenarios.

There are 54 unique genome IDs overall, not 600 independent samples. The accepted
count is 572 under declared comparator scopes; 24 clinical vep-rs timing-only
and four SnpEff strict failures are retained and excluded from validated claims.

## Figures

**Preferred final exports: `figures-v2/`.** The earlier `figures/` set remains
unchanged; v2 adds the payload-bound build scenario and expanded reusable data.

| Figure | Purpose | PDF |
|---|---|---|
| 1 | Matched assay speedups and absolute WGS savings | [Figure 1](figures-v2/manuscript/Figure1_assay_performance.pdf) |
| 2 | 52 WGS genomes × three cache strategies | [Figure 2](figures-v2/manuscript/Figure2_wgs_cache_strategies.pdf) |
| 3 | Measured stages and conditional optimization headroom | [Figure 3](figures-v2/manuscript/Figure3_optimization_headroom.pdf) |
| S1 | Configurable delayed pipelines, 80% and 90% hits | [Figure S1](figures-v2/supplements/FigureS1_pipeline_complexity.pdf) |
| S2 | All accepted assays by assembly | [Figure S2](figures-v2/supplements/FigureS2_all_assays_by_assembly.pdf) |
| S3 | Fresh paired CPU and aggregate memory | [Figure S3](figures-v2/supplements/FigureS3_cpu_and_memory.pdf) |
| S4 | Historical/current cached context, not causal A/B | [Figure S4](figures-v2/supplements/FigureS4_historical_cached_comparison.pdf) |
| S5 | Explicitly assumed build-cost model | [Figure S5](figures-v2/supplements/FigureS5_illustrative_amortization.pdf) |
| S6 | Payload-bound historical fastVEP build-cost scenario | [Figure S6](figures-v2/supplements/FigureS6_recorded_build_scenario.pdf) |
| D1 | Unvalidated clinical vep-rs timings, diagnostic only | [Diagnostic D1](figures-v2/diagnostics/DiagnosticD1_veprs_recipe_timings.pdf) |
| Repository | Simple WGS runtime-reduction bars | [Repository graphic](figures-v2/repository/Repository_WGS_time_reduction.pdf) |

Every graphic has vector PDF/SVG and 600-dpi PNG. Composite manuscript figures
and S1 also have individual panel exports under `manuscript/panels/` and
`supplements/panels/`. Figures S2–S6 use faceting within one plot, not separately
lettered compound panels. Captions, not embedded titles/footnotes, explain scope.

## Reproduce locally

From the repository root, using the existing Python environment with project
dependencies and R with ggplot2, patchwork, scales and tibble:

```bash
.venv/bin/python -m benchmarks.performance_sprint.audit_final \
  benchmarks/performance_sprint/final-analysis-20261029

# Must be a NEW figure directory: previously exported graphics are never replaced.
Rscript --vanilla benchmarks/performance_sprint/analyse_final.R \
  benchmarks/performance_sprint/final-analysis-20261029 figures-new-review

.venv/bin/python -m pytest tests/test_final_analysis.py -q
Rscript --vanilla tests/test_final_analysis.R
```

The R entry point adds `~/.local/lib/R/4.5` to its search path; override with
`VCFCACHE_R_LIBRARY` if needed. No package installation or network request occurs
during analysis. Summary tables/RDS are derived artifacts and are regenerated;
the original evidence archive and existing figure files remain unchanged.
The run environment is recorded in [R_sessionInfo.txt](source/R_sessionInfo.txt).

To explore without rerendering:

```r
data <- readRDS("benchmarks/performance_sprint/final-analysis-20261029/tables/manuscript_data.rds")
names(data)
data$primary_summary
data$measurements
data$comparator_audit
data$target_headroom
```

The read-only remote collector is [snapshot_final.py](../snapshot_final.py).
It refuses to replace a completed source snapshot, verifies manifest READY
hashes and terminal comparisons, and never launches annotation. Recreating a
remote snapshot requires a new output directory and access to the retained
controller campaign paths. Full BCFs are intentionally not part of this compact
Git archive.
