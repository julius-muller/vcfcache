# VCFcache: implementation decision workspace

This is a **flat folder**: one Excel workbook, PNG-only figures, two build
scripts, a verifier, and the data needed to inspect/recreate the figures. No
new annotation was run. Existing figures, benchmark sources and application
code were not changed.

## Open these first

1. **[decision_runtime_comparison.xlsx](decision_runtime_comparison.xlsx)** —
   98 explicit configurations, identical rows/order in **Panel**, **WES**, and
   **WGS**. The current candidates are first; historical tuning is retained
   further down. Each annotation recipe has runtime seconds, paired speedup,
   sample count and validation status. Use the **+** groups to reveal min/max,
   p95, direct time, seconds saved, CPU time, memory, counts and actual backends.
   Settings/provenance columns D–F are also expandable. Filters and frozen
   headings are enabled. Blank cells genuinely mean no measurement.
2. **`01_cohort_recipe_speedup_{Panel,WES,WGS}.png`** — the broadest comparable
   evidence: Python-adaptive versus Rust-adaptive on the common 12-input
   GRCh38 assay set. Individual inputs remain visible beneath the median lines.
3. **[05_python_vs_rust_same_conditions.png](05_python_vs_rust_same_conditions.png)** —
   the marginal benefit of Rust rather than the whole caching benefit. Above
   1 means Rust's cached run was faster. This is historical matched-condition
   evidence, **not** a newly randomized same-node A/B test.
4. **[06_rust_125_incremental_gain.png](06_rust_125_incremental_gain.png)** —
   separates external-tool, linked-library, and input-preservation changes.
   Each point uses the same block's **own engine before 1.25**, not adaptive
   versus resident relabelled as a library improvement.

## The other PNG views

| Prefix | View | Assay-specific files |
|---|---|---|
| `02_historical_branch_speedup_` | Earlier main, legacy-at-head, original fused, adaptive, resident and explicit adaptive fallback | Panel, WES, WGS |
| `03_observed_cost_speedup_` | Speedup against observed direct microseconds/variant on log axes, faceted by evidence set | Panel, WES, WGS |
| `07_latency_vs_direct_` | Absolute direct versus cached seconds; below the labelled parity diagonal is faster | Panel, WES, WGS |
| `04_all_configurations_overview` | Zoomable heatmap of all configurations with measured base recipes | Combined |

There are **15 PNGs**: four assay-specific figure families (three files each)
plus three overall comparisons. All are generated with **R/ggplot2**. Branch
labels match the workbook, including configuration/CPU keywords. The extensive
heatmap is a navigation aid; the workbook is the readable detailed view.

## Reading the numbers correctly

- **Runtime** is measured complete-CLI elapsed seconds, excluding cache
  build/download and explicitly untimed preparation/auditing. It is not a
  single internal stage. `Sum measured wall` is the sum of exported timings,
  **not** a forecast of cohort turnaround or wall time on a parallel cluster.
- **Speedup** is each measurement's own recorded direct time / cached time.
  Medians are medians of those ratios, not a ratio of median times. Ratios below
  1 mean caching was slower. Green/red is descriptive, not a significance test.
- Runtime colours compare configurations **within one evidence set and recipe**.
  Speedup cells are green above 1.05 and red below 1; the small neutral region
  does not establish an equivalence margin statistically.
- The three sheets have the **same configurations and recipe columns**, even
  for an entirely unmeasured assay. No interpolation fills the gaps. A 1,000-
  genome/gnomAD-overlap pilot and cache-derived all-hit microbenchmarks are not
  presented as real assay performance.
- **98 rows do not mean 98 branches or independent experiments.** Rows preserve
  branch + measured configuration + evidence set. One implementation can have
  CPU, worker, decoder, memory, region, cache or toolchain variants. Known
  experiments with no accepted assay runtime, such as failed whole-cache RAM
  loading, retain a labelled blank row rather than a loading-time “speedup”.
- `Samples` counts distinct genome IDs. `Sample-node units` also distinguishes
  repeated access to the same genome on different workers. Technical repeats
  are reduced to their median within sample/node before distribution summaries.
  Missing worker IDs remain missing; repeated same-input records then collapse
  conservatively into one unit. p95 is blank below five observed units. Min/max
  and p95 are descriptive, **not** confidence intervals or extra replication.
- Full-cohort, one/two-input pilot, CPU budget, assembly, cache strategy and
  profiling conditions are kept separate. Different recipes/assemblies can add
  different information. Lines connect measured recipe summaries, **not** a
  causal dose-response fit. Observed direct microseconds/variant includes
  startup, parsing, I/O and output; it is not pure annotator computation.
- Native requests that ran adaptive are **separate fallback rows/series**.
  They are not resident WGS results. Ordinary “native” and Rust-adaptive do not
  refer to the same implementation.
- Hollow points and `*` heatmap values mean known annotation differences or
  timing-only evidence, **not** statistical significance. Clinical vep-rs WES/
  WGS and historical VEP HGNC_ID differences remain visible. Accepted cohort
  VEP measurements retain their documented comparator scope, not a new claim
  of strict byte identity. Failed SnpEff audits remain flagged. The workbook
  preserves successful execution timings even when the annotation audit failed.
- Main + 1.25 includes timed identity/duplicate validation absent from the
  trusted-input pilots; many cross-campaign comparisons also use different
  workers/dates/cache states. Their gap is **not** a clean engine-only effect.
  sl-w5 is AMD in historical campaigns; recent toolchain pilots exclude it.
- Profile CPU-scaling exports without a matched-CPU direct denominator keep
  elapsed times but leave speedup blank. The compact sprint export similarly
  lacks direct elapsed times; no denominator is guessed from another job.
- Peak memory is the exported **cgroup/Slurm-step peak, including page cache**,
  not necessarily process RSS. Missing resources are not invented.

## Scope and source traceability

`configurations.tsv` lists every row and its precise scope. `measurements.tsv`
retains 1,664 exported timing observations with condition IDs, worker, counts,
audit status, resource values and reference locations. Analysis strata overlap:
these are **not 1,664 independent genomes**. `summary.tsv` is the exact numerical
source of the workbook/plots. `sources.tsv` binds the original pinned Git exports
and the read-only recent-pilot snapshot by SHA-256.

`source_recent_pilots.json` contains compact checkpoints, resource/audit file
hashes and terminal states from the main, Python and Rust 1.25 pilots. It does
not contain BCFs. `BUILD.json` binds the workbook; `VALIDATION.json` records the
checks. The generator uses pinned commits, not arbitrary dirty worktree data.

No single grand ranking is computed across unrelated recipes or cohorts. The
point is to decide which implementation is worth maintaining, not choose the
fastest historical observation per sample.

## RStudio and rebuilding

```r
decision <- readRDS("benchmarks/decision_20261109/decision_data.rds")
names(decision)
decision$measurements
decision$summaries
decision$historical_python_rust_pairs
decision$recent_rust_upgrade_pairs
```

From the repository root (Python dependency is isolated by uv, not installed
into the project environment):

```bash
uv run --no-project --with openpyxl==3.1.5 --python .venv/bin/python \
  benchmarks/decision_20261109/build_decision.py
Rscript --vanilla benchmarks/decision_20261109/figures_decision.R benchmarks/decision_20261109
uv run --no-project --with openpyxl==3.1.5 --python .venv/bin/python \
  benchmarks/decision_20261109/verify_decision.py
```

Default rebuild uses the local recent-pilot snapshot, requiring no SSH. Explicit
`--refresh` re-reads those immutable remote reports; it does not run benchmarks.
Only this decision folder's derived exports are replaced on rebuild. Previously
printed publication figures are never overwritten. R dependency/session details
are in `R_sessionInfo.txt`; all plotted data are saved as named tibbles in the RDS.
