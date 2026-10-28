# Performance sprint and manuscript rerun — deadline 2 November 2026

Target delivery: **Monday 2 November, 18:00 CET**. Performance work takes
priority in the first three days; candidate freeze is **no later than** Friday
30 October at noon. Freeze and submit earlier as soon as the checks pass; do not
leave ready work idle merely to follow a calendar date.
This file is the handoff entry point, not a claim that future stages ran.

## Feedback checkpoints and cluster utilization

**28 October header recovery:** user-approved fix deployed as `787a3081a` /
helper 0.7.1 after the adaptive guard rejected compatible receiving `AC Number=.`
headers. Real repaired gate 2001 passed strict complete-record/header comparison
on 165 variants with zero mismatches. Jobs **2004–2038** replace only the remaining
35 affected conditions; successful runs and all original data remain untouched.
Three completed WGS direct controls are reused (4.8 node-hours preserved), and
the fourth in-flight control is protected by an explicit predecessor dependency.
Final collector **2040** waits for both original collector 2000 and repair
collector 2039. All six workers remain active. The main checkpoint tree retains
all 600 logical slots, with repairs in a003/a004 and explicit build provenance.
Details and exact paths: `benchmarks/performance_sprint/ACTIVE_MANUSCRIPT_CAMPAIGN.md`.

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
- **User selected A on 27 October:** stop optimization work, freeze the supported
  candidate and execute the complete manuscript campaign. The 600-condition
  manifest is materialized: 528 unique historical cached conditions plus 72 new
  matched vep-rs pairs. Twelve identical historical WGS rows are aliases, not
  repeated measurements. Primary engine: adaptive; native remains opt-in.
- Preparation job 1376 exposed an independent C2 PGP sample-set lookup and failed
  before measurements. Corrected preparation 1377 succeeded. First primary jobs
  1378–1383 yielded three completed checks and three blocked checks; automatic
  dispatcher 1385 correctly refused to release the bulk campaign. No failure is
  silently discarded. Receipts preserve both immutable deployments.
- Two clinical checks exposed limitations in historical projection-only digests.
  Untimed overlay audit 1386 reproduced the new outputs with the legacy INFO
  overlay. Fresh direct controls 1387 confirmed full fastVEP equality; SnpEff
  additionally exposes reference-ID enrichment and one non-ASCII clinical-string
  encoding difference. These are validation/configuration issues, not evidence
  of a routing speed regression. Do not waive strict complete-record checks.
- Recovery keeps original caches untouched. Canonical transport copies repair
  drifted worker-local replicas into private node-local directories, outside
  timing, without changing the original local/shared storage class. New clinical
  pairs replace the 180 projection-only baselines with full retained direct
  outputs. The SnpEff clinical recipe explicitly preserves input IDs (`-noId`)
  and pins its historical ASCII text-encoding contract in both arms; the existing
  accented-string limitation must not be represented as lossless Unicode output.
  These required direct controls add approximately 53 node-hours (about nine
  hours across six workers), not routine technical replicates.
- Unchanged, already-complete first checks can be carried forward only with
  matching product files, inputs, cache bytes, parameters, CPU/engine settings
  and baseline binding. Changed configurations and failed checks need explicit
  fresh validation. Known clinical vep-rs differences remain timing-only and
  retain both full outputs; no upstream tool changes are authorized.

**Live launch verified at 19:33 UTC / 20:33 CET on 27 October:** preparation
1398 completed; all six primary conditions passed (two unchanged successful
checks reused). Dispatcher 1405 released **594 jobs, 1406–1999**, with afterany
collector **2000**. All six workers were running, with 588 annotation jobs queued
and zero failures recorded. The full 600-condition campaign is now executing;
this is not a claim that it has completed. Host and Docker suites both pass
**646 tests, 21 skipped**; GitHub CI passed deployed revision `fe197f34d`.

No new user decision is currently needed. Next operational check-in:
**22:35 CET on 27 October**, approximately two hours after bulk release. Slurm
and final collection advance automatically; unexpected failures still require
active inspection. Provisional completion: **30–48 hours after bulk release**,
to be revised from measured throughput. See
`benchmarks/performance_sprint/ACTIVE_MANUSCRIPT_CAMPAIGN.md` for exact paths and
days-later review commands. Do not resume closed optimization experiments.

Commands and launch receipts live in `benchmarks/performance_sprint/`.
The live Slurm campaign, not this file's age, determines current progress.
There is no implied background agent: Slurm jobs and collectors persist without
the conversation, but code-review/freeze decisions require a subsequent work turn.
