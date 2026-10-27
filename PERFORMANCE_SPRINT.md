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

## Current implementation

Branch: `feature/performance-sprint-20261027`; base `fd805f0b3`.

- Optional helper 0.5.0 adds an explicit indexed-miss mode. Standard behavior
  remains plain BCF unless the experiment is requested.
- Environment-only screening switches:
  `VCFCACHE_EXPERIMENT_INDEXED_MISSES=0|1` and
  `VCFCACHE_EXPERIMENT_WINDOW_BP=1024|4096|16384`.
  Both are frozen per cell and reported in statistics/command metadata.
- The first screening matrix has **14 cached-only cells** on A06/C05: indexed
  spooling for core WGS and fastVEP WES, plus native WES windows and controls.
  Direct and uninstrumented cached references come from the verified original
  campaign; SHA256 checks precede timing. No direct runs are scheduled by it.
- Changes are experimental and **not yet accepted by performance gates**.

Commands and launch receipts live in `benchmarks/performance_sprint/`.
The live Slurm campaign, not this file's age, determines current progress.
There is no implied background agent: Slurm jobs and collectors persist without
the conversation, but code-review/freeze decisions require a subsequent work turn.
