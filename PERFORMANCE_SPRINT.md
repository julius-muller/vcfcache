# Performance sprint and manuscript rerun — deadline 2 November 2026

Target delivery: **Monday 2 November, 18:00 CET**. Performance work takes
priority in the first three days; candidate freeze is **no later than** Friday
30 October at noon. Freeze and submit earlier as soon as the checks pass; do not
leave ready work idle merely to follow a calendar date.
This file is the handoff entry point, not a claim that future stages ran.

## Feedback checkpoints and cluster utilization

- Every status/handoff must state: **what is running; what advances automatically;
  what requires active agent work; whether user feedback is needed; and the next
  approximate check-in/resumption time**. Distinguish elapsed job ETA from the
  additional active work needed for preparation. Say explicitly if nothing is
  queued and the plan will not advance after the turn ends.
- Do not ask for renewed permission for already-approved technical work. Request
  user feedback only for a genuine choice, scope change or missing authority;
  explain the decision and when it is needed. If a later active turn is required,
  call it a **resumption point**, not a request for scientific feedback.
- During benchmark execution, queue all ready independent conditions across all
  healthy eligible workers, with enough queued work to refill freed nodes without
  another conversation turn. Submit longest conditions first with realistic wall
  limits; use dependencies for already-known subsequent stages where safe.
- Each idle or excluded worker needs a concrete reason in status reports:
  preparation/validation not finished (with remaining work), fewer remaining
  independent cells than nodes, resource/health failure, or a required controlled
  measurement. Re-check exclusions instead of inheriting them indefinitely.
- Keep matched per-condition CPU/memory settings and isolated timed steps. Do not
  change 8-CPU conditions to 16 CPUs, pack competing timing cells on one node, or
  create unnecessary repeats merely to improve an apparent utilization number.
  Slurm's current benchmark partition is node-exclusive.
- `sl-w5` passed a fresh execution/affinity/resource preflight on 27 October.
  The manuscript runner no longer hard-codes its exclusion. Any future exclusion
  must be frozen with a named node and a written reason. Recheck health before
  submission; use all six workers when they remain eligible.

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
3. **By Friday noon, earlier if ready:** freeze selected changes, packaging, defaults, helper and
   recipes after host/container/strict-record testing. Native remains opt-in.
   Do not describe an experimental switch as an accepted product improvement.
4. **Immediately after freeze, targeting completion by Sunday:** rerun cached manuscript cells using approved direct
   references; separate jobs for long direct exceptions. One Slurm step/cgroup
   per condition, explicit retries, durable per-cell state, afterany collectors.
   Target all measurements/audits by Sunday 18:00. Monday: final consolidation,
   R/ggplot2 figures in a new directory and handoff, not silent sample reductions.

## Current implementation (27 October — engineering screen complete)

Branch: `feature/performance-sprint-20261027`; base `fd805f0b3`.

- Optional helper **0.7.0** supports indexed miss spools, native query windows,
  trusted-hit assembly and a bounded native router inside the existing adaptive
  workflow. Python adds a bounded decoded-payload memory preflight. Every switch
  is experimental, frozen per condition and disabled by default.
- **36/36 engineering cells completed**, all strict direct/frozen-cache output
  checks passed: initial 1347–1350; trusted-hit
  1353–1356; memory preflight 1359–1360; bounded route 1362–1365. Each phase has an
  afterany collector. Final reversed-order memory-preflight confirmation is
  **job 1368**, collector **1369**, campaign
  `/results/performance-sprint/20261027T112225Z-2ea9a8cc6`.
  The engineering budget is exhausted: no extra exploratory matrix is authorized.
- Initial and trusted-hit screens completed with strict direct and frozen-cache
  parity. Neither indexed spooling, smaller windows nor trusted-hit assembly
  showed a consistent >=5% gain across the two genomes: **do not promote them**.
- Memory preflight: A06 **62.68 → 45.90 s (26.8% reduction)**, C05
  **79.32 → 62.22 s (21.6%)**, all strict comparisons pass. Both controls otherwise
  allocate an 8-GiB resident buffer and fall back. The prediction samples at most
  128 records and takes <0.4 s here; its size estimate is heuristic, not an RSS
  guarantee. Reversed-order A06 confirmation also passed: **65.33 → 45.60 s
  (30.2%)**. This supports candidate selection for avoiding failed resident loads,
  not a universal speedup over already-adaptive routing. Defaults remain unchanged
  until the release-candidate gate.
- Bounded routing: fastVEP WES **51.91 → 47.81 s (7.9%)** and
  **52.51 → 49.67 s (5.4%)**, strict parity passes. Core WGS is **23.3% / 31.0%
  slower** on A06/C05. **Do not make it a universal default**.
- Latest host suite: **633 passed, 21 skipped**. A freshly built standalone
  Docker test image also passes **632 tests, 21 skipped**, without mounting the
  source checkout; the small manuscript fixture tables are included explicitly.
  Ruff, Black and Rust formatting pass. The replay harness has a real-CLI
  integration test with only Slurm and the resource sampler mocked; this is not
  a claim of cluster replay execution.

## Manuscript rerun preparation

- `inventory-audited-20261027/` preserves **540 historical cached rows**, including
  all SnpEff conditions, 52-genome strategies and both synthetic-load cache sets.
  This is before provenance-based deduplication; no sample reduction is approved.
  Historical cached time sums to **108.2 node-hours** (~18 h over six workers),
  before audits/queueing and the pinned vep-rs extension. It is not a completion ETA.
- Read-only audits: job 1352 verified **23,200 annotation assets**, no failures.
  Job 1358 verified **64 input/recipe/region files**, no mismatches; its one missing
  unpadded WES interval file is not used by the manuscript. All 52 original genome
  hashes match. The fastVEP and two SnpEff wrappers also match their frozen copies.
- **Job 1367 finished**: 192 retained-output rows verified, 180 C2 digest-only
  rows, 168 unresolved reference checks. These are audit limitations, not failed
  benchmark runs. The 156 fastVEP lookups pointed to an obsolete campaign;
  `references.py` now derives the campaign from the SHA-bound accepted source
  `SNAPSHOT.json`. The 12 synthetic-load conditions retain timings and semantic
  metadata, but not per-delay direct BCFs; they are explicitly classified as
  metadata-only, not falsely promoted to complete-record output evidence.
- **Read-only supplemental job 1370 finished successfully**: 156 fastVEP
  retained-output rows verified and 12 synthetic metadata-only rows verified;
  **zero unresolved rows remain**. Across both audit passes, all 540 inventory
  rows now have explicit evidence classes: 348 retained-output, 180 C2 digest-only
  and 12 synthetic metadata-only. This is evidence classification, not automatic
  approval of timing reuse. No annotation was repeated. Request:
  `reference-request-v2-20261027.json`; receipt: `audit-supplement-20261027.json`.
  The synthetic conditions still require an explicit composed-reference
  decision before reuse; C2 digest-only limitations likewise remain visible.
- `manuscript.py` implements a cached-only replay runner with a candidate-freeze
  gate, exact task-list binding, per-task approved reuse evidence, matched CPU
  settings, cache/input/resource fingerprints, one independent Slurm step per
  cell, explicit retries and durable afterany reviews. Successful cells cannot
  be silently repeated. Killed/stale running states are not treated as success.
- **The complete replay task document, baseline approval and candidate freeze are
  not yet finalized; the manuscript rerun has NOT been submitted.** Pinned vep-rs
  core/clinical still needs its 12-genome matched-assay extension materialized.
  Preserve its known ClinVar discrepancies; no upstream tool work is authorized.

Next active work (no new user decision currently required): select only supported
changes, complete per-row reuse bindings (including CPU/preprocessing
and digest-only limitations), materialize the replay tasks and vep-rs extension,
then package/test/freeze before the Friday submission gate. If audits reveal a
baseline mismatch, rerun only that direct exception in a separate job.

At the 27 October afternoon handoff, no long-running job remains. The critical
path is active preparation, not Slurm. Estimate **1–3 hours of further active
work** to reach the next candidate/replay-readiness checkpoint; this is not an
unattended wall-clock promise or an automatic scheduled launch. Resume promptly
rather than waiting for Friday. Report a revised estimate if preparation exposes
additional binding/validation work.

Commands and launch receipts live in `benchmarks/performance_sprint/`.
The live Slurm campaign, not this file's age, determines current progress.
There is no implied background agent: Slurm jobs and collectors persist without
the conversation, but code-review/freeze decisions require a subsequent work turn.
