# Decision brief

29 October 2026. All approved runs are finished; no annotation job is queued.
The full evidence and limitations are in [ANALYSIS.md](ANALYSIS.md).

## Recommended decision: freeze and use the current dataset

The common GRCh38 cohort supports median WGS acceleration of **10.23× VEP,
10.81× clinical SnpEff, 2.88× clinical fastVEP and 1.32× core vep-rs**. Native
is worth retaining as opt-in, but the full campaign's primary is adaptive.
FastVEP and core-vep-rs WES losses remain visible; rich vep-rs WES/WGS ratios
remain timing-only. This is a strong workload-specific product story, not an
“always faster” claim.

## Decisions requiring user direction

| Priority | Choice | Recommendation | What it unlocks |
|---|---|---|---|
| 1 | Freeze dataset and move to manuscript integration? | **Yes** | Revise figures/methods around the complete campaign, without another whole rerun |
| 2 | Change release engine selection? | Review availability-based adaptive selection; retain native opt-in and legacy fallback | Ordinary users can obtain the measured path; current CLI default is still legacy |
| 3 | More engineering before publication? | **Optional, bounded WGS routing/assembly pilot only**; do not block the paper | Quantify an actual improvement before authorizing another full replay |
| 4 | Is a validated rich-vep-rs headline indispensable? | Not a release blocker; otherwise a separately scoped recipe-equivalence follow-up | Pinned clinical data cannot honestly substantiate that claim |

No decision is needed to let existing jobs finish: **they have already finished**.
No autonomous benchmark or release work will advance after this analysis turn.
The next user feedback is needed **now, before** new implementation, submissions,
default changes or manuscript integration; these are different scopes from the
approved consolidation. Nothing requires keeping the cluster busy merely to
avoid idle nodes.

## If priority 3 is chosen

1. Freeze one ordinary and one low-hit SGDP genome, core vep-rs and clinical
   fastVEP, the current output format and CPU/memory settings. Reuse compatible
   direct baselines, but run a small matched cache-side control/candidate check.
2. Inspect routing and assembly separately. For core-vep-rs WGS they cost about
   **20/18 s**, respectively; generic CLI overhead is only about **2 s**.
3. Require >=10% repeatable end-to-end gain, strict full-record/header parity,
   no >5% repeatable regression on the other workload, and an explicit memory
   cap. A stage-only microbenchmark is insufficient.
4. Stop after the agreed short pilot if the gate fails. If it passes, decide
   whether the change warrants targeted replay or a complete campaign. Do not
   start a full rerun while the engine decision is still open.

## What not to spend time on now

- Do not tune core-vep-rs WES with an expectation of a quick, spectacular win:
  parity needs roughly **77%** less cache-side time with miss annotation fixed.
- Do not promote small windows/indexed spools/trusted assembly: their screening
  results were inconsistent or regressive. Do not claim the native fallback
  guard improves already adaptive runs.
- Do not discard native's real WES benefit; keep its maintenance surface bounded
  through shared orchestration, explicit opt-in and parity/backend tests.
- Do not patch/upgrade the pinned vep-rs showcase, ignore failing fields, invent
  setup costs, or claim another full-cohort native evaluation was performed.
- Do not initiate a new format, database service or whole-engine rewrite without
  a demonstrated end-to-end bottleneck-specific benefit.

## Medium-term opportunity

Cohort miss-union/batching can reduce repeated costly work and invocations, with
an equally batched uncached control. It is more promising for heavy workflows
than shaving their already small wrapper fraction. Safe subset-invariant recipes,
per-sample genotype/metadata preservation, merge cost and versioned provenance
are prerequisites—not optional follow-up checks.
