# Decision figures v2: historical Intel-only sensitivity analysis

**No benchmarks were rerun or submitted.** This is a reanalysis of existing
measurements. The entire [v1 folder](../decision_20261109/00_START_HERE.md) is
unchanged, including its figures and Excel workbook.

## Open these first

- [Excel comparison](decision_runtime_comparison.xlsx): Panel, WES and WGS,
  with the same **98 configuration rows and 14 recipe blocks** as v1.
  Expand the column groups for runtime min/max/p95, direct runtime, sample
  counts, resources and validation status. Blank now also means excluded.
- [Figure 5: Python versus Rust](05_python_vs_rust_same_conditions.png): the
  requested comparison after exclusion. Colour distinguishes the same Intel
  worker from different Intel workers; both are historical comparisons, not
  randomized interleaved A/B trials.
- [WGS speedups](01_cohort_recipe_speedup_WGS.png),
  [WES speedups](01_cohort_recipe_speedup_WES.png), and
  [Panel speedups](01_cohort_recipe_speedup_Panel.png): medians and individual
  inputs, with the retained **n** shown for every recipe.

All **15 PNGs** were regenerated with R/ggplot2. The other figure families
retain the v1 names: historical branch matrix (`02`), observed direct cost
versus speedup (`03`), all-configuration overview (`04`), Rust toolchain
incremental gain (`06`), and direct versus cached latency (`07`). The folder
is flat; it contains no subfolders, PDFs or SVGs.

## What was excluded

1. A timing is excluded if its cached run used historical **sl-w5**.
2. A ratio is excluded if its **direct denominator** used historical sl-w5,
   even when the cached run was on Intel.
3. Placement must be established from pinned worker/resource metadata or a
   duration-matched direct-reference resource record. Unknown placement is
   excluded, not presumed Intel. Logical profile groups such as `g00` are
   resolved to actual workers.
4. Paired Python/Rust cohort views use the **intersection** of eligible inputs:
   if either engine or the direct denominator is excluded, both views are
   removed. This avoids comparing different sample compositions between
   engines within a recipe/assay. Different recipes can still retain different
   sample sets; their medians are not controlled cross-tool comparisons.

Of 1,664 exported analysis observations, **876 remain and 788 are excluded**:
552 explicitly involve sl-w5; another 73 have unresolved placement without an
explicit sl-w5 assignment; 163 otherwise eligible partner observations are
removed to keep the cohort pairing balanced. These are overlapping analysis
views, **not** counts of independent genomes or new runs.

The main cohort retains **88 matched Python/Rust pairs**, with **4–9 samples
per recipe/assay**. Every configuration and recipe column remains in Excel,
including completely empty configurations. In particular, the synthetic-delay
direct times are composed references whose component worker provenance is not
bound by the exported baseline path. They are excluded here rather than
assigned the ordinary VEP reference's worker. Older resident/projection and
engineering-sprint cells without recoverable worker evidence are also blank.
Their original measurements remain in v1 and the exclusion export.

No runtime is rescaled, corrected, substituted with a different direct arm,
or selected for being fast. Technical repeats are reduced within sample/worker,
as in v1. Min/max and p95 are descriptive, not confidence intervals; p95 is
blank below five sample/worker units. Annotation-audit caveats remain visible;
hardware filtering does not convert timing-only evidence into validated output.

## The suspicious SnpEff WGS split

The original approximately 0.5× / 1× / 2× clusters in Figure 5 largely track
the hardware pairing. After removing sl-w5 from both cached arms and direct
references, **four SnpEff clinical WGS pairs remain**, with Python/Rust cached
time ratios approximately **1.02–1.09×** (median **1.049×**). Above 1 means
Rust was faster. The twofold engine advantage/disadvantage is not reproduced
in this Intel-only subset.

This supports hardware confounding; it does not identify CPU, storage,
filesystem caching or shared-host contention as the sole cause. It also does
not establish equal performance among the remaining Intel workers.

## Relation to the 9 November Slurm audit

The [audit snapshot](SLURM_CPU_PLACEMENT_AUDIT_20261109.md) was read from
`ubuntu@10.133.242.60:/home/ubuntu/portal-services/SLURM_CPU_PLACEMENT_AUDIT_20261109.md`.
Its final updates supersede its initial inspection/recovery recommendations:

- Historical sl-w5 exposed AMD EPYC-Rome; the other workers exposed Intel
  Broadwell. The audit records sl-w5 being **replaced by an Intel VM**, retaining
  its data volume.
- The audit then records **all seven nodes upgraded/rebooted**, with kernel
  7.0.0-43-generic and Apptainer 1.5.5. The historical timings predate this
  environment; the new sl-w5 identity must not relabel the old measurements.
- A repeated SHA-256 microbenchmark varied by up to about **15% single-threaded
  and 8% at 16 threads on the same worker**. Those are audit microbenchmark
  observations, not uncertainty estimates for these annotation timings.
- Matching guest CPU labels do not establish identical underlying host
  performance. CPU placement remains unconstrained for future rebuilds.

Consequently this is a **sensitivity analysis**, not a repaired publication
campaign or causal engine ranking. Excluding the conspicuous mixed-hardware
worker is useful now; rerunning is a separate subsequent decision.

If we later investigate residual differences, start with a small, paired,
interleaved same-node A/B check of the disputed WES/WGS recipes, repeat it
enough to quantify timing spread, and record VM identity, CPU model, kernel,
container/tool versions, affinity, background load and storage placement.
Use both direct and cached controls in that new environment. Do not silently
combine new cached measurements with old pre-upgrade direct baselines.
**No such follow-up was launched in this step.**

## Data and reproduction

- `measurements.tsv`, `summary.tsv`, `configurations.tsv`: exact retained
  numerical inputs; `retained_counts.tsv`: before/after counts per cell.
- `excluded_measurements.tsv`: every excluded original observation, with
  resolved cached/direct workers, proof locations and the exclusion reason.
- `node_provenance.json`: frozen worker proofs, pinned metadata hashes, and
  compact read-only direct-reference metadata retrieved from the controller.
- `BUILD.json`: workbook hash, counts, policy, and hashes of every v1 file
  checked for preservation. `VALIDATION.json`: independent verification.
- `decision_data.rds`: named tibbles, including the filtered historical pairs.
  `R_sessionInfo.txt` records plotting dependencies.

From the repository root:

```bash
uv run --no-project --with openpyxl==3.1.5 --python .venv/bin/python \
  benchmarks/decision_20261109_v2/build_v2.py
Rscript --vanilla benchmarks/decision_20261109_v2/figures_decision.R benchmarks/decision_20261109_v2
uv run --no-project --with openpyxl==3.1.5 --python .venv/bin/python \
  benchmarks/decision_20261109_v2/verify_v2.py
```

Default rebuilding is offline and reuses v1's pinned evidence and summarization
code. `build_v2.py --refresh` reads historical controller JSON metadata over
SSH; it does **not** schedule or execute any benchmark. It is unnecessary for
reproducing these frozen figures. No Python package is added to the project
environment. Plotting requires the same R packages as v1.
