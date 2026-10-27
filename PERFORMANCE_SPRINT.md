# Performance sprint and manuscript rerun — deadline 2 November 2026

Target delivery: **Monday 2 November, 18:00 CET**. Performance work takes
priority in the first three days; candidate freeze is Friday 30 October at noon.
This file is the handoff entry point, not a claim that future stages ran.

## Fixed decisions

- Pinned vep-rs core/clinical is a speed showcase; do not patch or upgrade it.
  Preserve inherited clinical discrepancies and never claim exact equivalence.
  New cache-engine differences against the frozen cached reference still fail.
- 2x core-vep-rs WGS is a stretch target, not a release gate.
- Reuse uncached timings and retained outputs when provenance/settings match.
  Missing historical full outputs are an audit limitation, not permission to
  claim a full-record comparison: distinguish accepted historical digests from
  new engine-versus-legacy full-record parity checks.
- Preserve all original manuscript experiment families and biological samples,
  including SnpEff, cohort caches and synthetic loads. No routine technical repeats.
- No whole-engine rewrite, new storage format, annotator integration, or public
  package release during this sprint. Default selection stays unchanged until
  the measured acceptance/release gate.

## Execution order and gates

1. **Tuesday:** inventory the complete manuscript matrix, bind historical direct
   references to inputs/recipes/resources/CPU settings, and record reuse decisions.
   Do not rerun a direct condition merely because a cached condition failed.
2. **Tuesday–Thursday:** screen indexed miss spooling and 1/4/16-kb native
   windows; investigate early fallback only after these measurements. At most
   36 cached timed controls, including confirmations and regression checks.
   Adopt only reproducible >=5% end-to-end improvements without >5% reproducible
   regressions, or removal of substantial demonstrated fallback waste.
3. **Friday noon:** freeze selected changes, packaging, defaults, helper and
   recipes after host/container/strict-record testing. Native remains opt-in.
   Do not describe an experimental switch as an accepted product improvement.
4. **Friday–Sunday:** rerun cached manuscript cells using approved direct
   references; separate jobs for long direct exceptions. One Slurm step/cgroup
   per condition, explicit retries, durable per-cell state, afterany collectors.
   Target all measurements/audits by Sunday 18:00. Monday: final consolidation,
   R/ggplot2 figures in a new directory and handoff, not silent sample reductions.

## Current implementation (27 October, 12:20 CET)

Branch: `feature/performance-sprint-20261027`; base `fd805f0b3`.

- Optional helper **0.7.0** supports indexed miss spools, native query windows,
  trusted-hit assembly and a bounded native router inside the existing adaptive
  workflow. Python adds a bounded decoded-payload memory preflight. Every switch
  is experimental, frozen per condition and disabled by default.
- **34/36 engineering cells submitted**: initial 1347–1350; trusted-hit
  1353–1356; memory preflight 1359–1360; bounded route 1362–1365. Each phase has an
  afterany collector. The remaining two cells are reserved for reversed-order
  memory-preflight confirmation, not another exploratory matrix.
- Initial and trusted-hit screens completed with strict direct and frozen-cache
  parity. Neither indexed spooling, smaller windows nor trusted-hit assembly
  showed a consistent >=5% gain across the two genomes: **do not promote them**.
- Memory preflight: A06 **62.68 → 45.90 s (26.8% reduction)**, C05
  **79.32 → 62.22 s (21.6%)**, all strict comparisons pass. Both controls otherwise
  allocate an 8-GiB resident buffer and fall back. The prediction samples at most
  128 records and takes <0.4 s here; its size estimate is heuristic, not an RSS
  guarantee. Confirmation precedes adoption.
- Bounded routing: fastVEP WES **51.91 → 47.81 s (7.9%)** and
  **52.51 → 49.67 s (5.4%)**, strict parity passes. Core-WGS audits remain in
  progress; the first cell's CLI time is slower. **Not a universal default**.
- Host/container suite before the new replay harness: **618 passed, 21 skipped**
  in each. The replay harness has a real-CLI integration test with only Slurm and
  the resource sampler mocked; this is not a claim of cluster replay execution.

## Manuscript rerun preparation

- `inventory-audited-20261027/` preserves **540 historical cached rows**, including
  all SnpEff conditions, 52-genome strategies and both synthetic-load cache sets.
  This is before provenance-based deduplication; no sample reduction is approved.
  Historical cached time sums to **108.2 node-hours** (~22 h over five workers),
  before audits/queueing and the pinned vep-rs extension. It is not a completion ETA.
- Read-only audits: job 1352 verified **23,200 annotation assets**, no failures.
  Job 1358 verified **64 input/recipe/region files**, no mismatches; its one missing
  unpadded WES interval file is not used by the manuscript. All 52 original genome
  hashes match. The fastVEP and two SnpEff wrappers also match their frozen copies.
- **Job 1367** hashes retained complete outputs against historical metrics,
  checkpoints every row and explicitly identifies the **180 C2 rows with only
  annotation/key digests retained**. It does not rerun annotation or automatically
  approve timing reuse. `reference-request-20261027.json` resolves all 540 rows.
- `manuscript.py` implements a cached-only replay runner with a candidate-freeze
  gate, exact task-list binding, per-task approved reuse evidence, matched CPU
  settings, cache/input/resource fingerprints, one independent Slurm step per
  cell, explicit retries and durable afterany reviews. Successful cells cannot
  be silently repeated. Killed/stale running states are not treated as success.
- **The complete replay task document, baseline approval and candidate freeze are
  not yet finalized; the manuscript rerun has NOT been submitted.** Pinned vep-rs
  core/clinical still needs its 12-genome matched-assay extension materialized.
  Preserve its known ClinVar discrepancies; no upstream tool work is authorized.

Next review: finish the core-WGS comparisons and final confirmation; select only
supported changes, complete per-row reuse bindings (including CPU/preprocessing
and digest-only limitations), materialize the replay tasks and vep-rs extension,
then package/test/freeze before the Friday submission gate. If audits reveal a
baseline mismatch, rerun only that direct exception in a separate job.

Commands and launch receipts live in `benchmarks/performance_sprint/`.
The live Slurm campaign, not this file's age, determines current progress.
There is no implied background agent: Slurm jobs and collectors persist without
the conversation, but code-review/freeze decisions require a subsequent work turn.
