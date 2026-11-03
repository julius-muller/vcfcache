# Rust WES pilot — 3 November 2026

## Decision

**Partial success, not an all-sample solution.** Genuine cached fastVEP clinical
WES speedups are achievable without retaining the whole cache in RAM. Generic
worker allocation also brings vep-rs core to approximate median break-even,
but does not make caching faster on every core-WES input. Neither target shows
a consistent 2× gain. No direct-annotation bypass, product-default change,
assay-specific dispatch, or new publication campaign was introduced.

The best-supported fastVEP choice here is the **existing resident-input path
with four lookup workers**, not the new CSI-coalescing prototype. Do not present
this pilot as proof that the new prototype itself delivered the whole gain.

## Cold application-start results

Speedup is each input's measured direct CLI time divided by its cached CLI time.
All configurations have the same eight-CPU affinity allocation. Each recipe has
12 independent WES inputs, 74,215–108,539 records, and cache hits of 71.1–95.5%.
The inputs, recipes, and caches are taken unchanged from the frozen manuscript
matrix; no easier subset was substituted.

| Cached configuration | fastVEP clinical: median speedup | Faster inputs | vep-rs core: median speedup | Faster inputs |
| --- | ---: | ---: | ---: | ---: |
| Rust-adaptive control | 0.92× | 4/12 | 0.57× | 0/12 |
| Resident input, four lookup workers | **1.24×** | **12/12** | 0.80× | 2/12 |
| CSI-coalesced resident input, four workers | 1.25× | 12/12 | 0.78× | 3/12 |
| Resident input, eight lookup workers | 1.25× | 10/12 | **1.04×** | **7/12** |

Four-worker resident-input fastVEP speedups range from 1.11× to 1.92×. Every input
clears the prespecified 5% margin. Coalescing has a similar median but is faster
than ordinary resident input in only 5/12 fastVEP cases; it is not a reliable
incremental win. It reduces query count slightly, not the dominant recipe cost.

Eight workers improve cached vep-rs core time over four workers in 10/12 inputs,
with a median paired gain of 1.17×. Core lookup time falls from a median 3.34 s
to 1.96 s. Core speedups against direct span 0.56–2.14×, however: five inputs
still lose. Eight workers regress seven fastVEP cached times relative to four
workers, and two fastVEP inputs become slower than direct. There is **no single
configuration here that achieves all-sample break-even for both recipes**.

These are one-run-per-input/configuration engineering observations, not
publication-quality effect estimates or evidence of linear parallel scaling.
The eight-worker configuration was selected on two difficult core-WES cases,
then tested on the remaining ten core inputs and all twelve fastVEP inputs.
Those first two timings are retained, not repeated or replaced by faster runs.

## What remains expensive

With eight lookup workers, the measured miss-recipe invocation occupies a
median 69% of core cached elapsed time and 82% of fastVEP cached elapsed time.
This stage includes the annotator process, database loading, recipe commands,
and recipe-owned I/O; it is **not** a measurement of pure per-miss annotation.
Median miss-recipe elapsed times are 11.44 s for core and 36.13 s for fastVEP.

Eliminating all other measured work while holding those invocation times fixed
would give median ratios of 1.72× and 1.58×, respectively. These are conditional
optimistic counterfactuals, **not universal ceilings**: recipe timing varies
between invocations, and this experiment does not estimate a fixed startup
constant or establish that miss annotation is intrinsically slower than direct.

There is a plausible generic decoder interaction worth isolating later. At
eight CPUs, four lookup workers share four extra BGZF codec threads; eight
lookup workers leave no extra codec pool. The bundled HTSlib enables BGZF_CACHE,
but its multithreaded BGZF read path bypasses the ordinary single-threaded block
cache lookup. Thus this screen changes both query parallelism and decoder mode.
It cannot attribute the gain to worker count alone. A future discriminating
experiment should hold lookup workers fixed while changing codec-pool use,
then examine query-window granularity. This can remain region-driven, not
assay- or annotator-specific. That experiment was not added to this bounded run.

Panel and WGS successes do not imply that every intermediate-size input must
win: sparse lookup, startup, compression, and saved annotation work scale
differently. The observed core-WES result is not evidence of an output defect.

## Whole-cache RAM residency: stopped on capacity evidence

A separate prototype loads immutable standard-BCF records once, then uses
per-contig position search and exact allele matching across samples. Unit tests
cover two-job reuse with zero, partial, and full hit configurations. The first
implementation eagerly unpacked INFO. A second retained lazy BCF records and
unpacked only matching clones. Neither representation fits the pilot budget.

| Recipe | Representation | Records before stop | Loading time | Largest helper RSS |
| --- | --- | ---: | ---: | ---: |
| vep-rs core | Eager | 9,502,720 | 31.53 s | 20.01 GiB |
| fastVEP clinical | Eager | 6,373,376 | 33.99 s | 20.05 GiB |
| vep-rs core | Lazy BCF | 10,305,536 | 27.44 s | 20.02 GiB |
| fastVEP clinical | Lazy BCF | 6,701,056 | 25.16 s | 20.04 GiB |

The helper samples its RSS every 16,384 loaded records, so small overshoot is
expected. The 20-GiB guard is **not** an aggregate hard limit; each Slurm job
has a separate 28,000-MB allocation. Target-record buffers are capped at 2 GiB.
The two twelve-input session attempts failed during cache loading. Their direct
controls and failure evidence are retained, but there are **no accepted warm
latency or amortization results**. This rejects these representations on these
workers, not every possible RAM representation or larger hardware.

## Validity, failures, and scope

- 124 accepted timed CLI conditions: 24 fresh direct baselines and 100 cached
  conditions/configurations. Every accepted cached output passes complete-record
  and relevant-header comparison; zero fields are ignored and no native fallback
  is accepted. Fresh direct outputs also match their retained frozen references.
- Timing includes the complete annotation CLI, normalization, lookup, miss
  recipe, output construction/compression, and indexing. Input extraction and
  semantic comparison are outside timing; pre-existing cache build costs are
  not included in single-sample latency.
- Separate Slurm steps make condition-level cgroup peaks meaningful. Memory
  counters include page cache; they are not silently relabelled as helper RSS.
- The deterministic screen has three inputs per recipe, initially direct-first.
  The other nine use counterbalanced condition order. Later configurations reuse
  these exact measured direct baselines and compare outputs against SHA-bound
  retained direct files. They do not rerun direct annotation for historical
  symmetry. Processes start fresh, but filesystem caches are not flushed.
- Direct and recipe timings vary substantially; especially for fast annotators,
  these ratios should not be interpreted as noiseless hardware constants. No
  sample-specific best-of-configuration result is reported.
- Six early coalesced attempts passed output comparison but failed their backend
  metadata audit because the controller omitted the coalescing marker. They were
  retained as failed attempts; after fixing the marker, six new audited attempts
  were accepted. No fastest-run selection was used.
- The recorded strict Bash script is now actually executed with `pipefail`.
  Count/index checks independently reject truncated outputs in the pilot.
- No VEP, other assays, annotations, new cache formats, automatic direct bypass,
  manuscript figures, or public defaults were changed. Previously printed figures
  remain untouched. All pilot jobs, 2673–2742, are terminal.
- Verification: **702 tests passed, 21 skipped**; targeted pilot tests, Ruff,
  Rust formatting, and Git whitespace checks passed.

## Recommendation

1. Retain the resident-input result as evidence that fastVEP clinical WES can
   genuinely benefit from caching. Do not promote CSI coalescing or eight workers
   universally on these results.
2. If continuing optimization, isolate codec-pool/block-cache behaviour before
   adding another lookup store or RAM representation. It is a smaller,
   annotator-independent change with an observed lookup-time signal.
3. For core-WES single-sample latency, also profile the miss recipe's loading/I/O
   separately. Lookup-only work cannot remove the large measured recipe fraction.
   Batch miss annotation is a possible **throughput** follow-up, not a substitute
   for single-sample speedup and not implemented as a shortcut here.
4. Keep the branch experimental. The focused pilot is finished; a full benchmark
   rerun or a product-default change requires a separate decision.

## Files and reproducibility

- `tables/cli_conditions.tsv`: 96 primary condition rows.
- `tables/worker_tuning.tsv`: 28 additional configuration rows, including one,
  two, and eight lookup workers; the eight-worker rows cover all 24 targets.
- `tables/cli_summary.tsv`, `worker_summary.tsv`, `capacity_probes.tsv`,
  `attempts.tsv`, and `pilot_data.rds`: analysis-ready data.
- `figures-v3/`: final R/ggplot2 plots; `figures-v2/` is retained unchanged.
- `terminal-evidence.tar.gz`, `EVIDENCE_INVENTORY.json`, and `SNAPSHOT.sha256`:
  immutable terminal metadata, manifests/READY hashes, commands, resources,
  comparisons, failed attempts, job IDs, and file integrity checks.
- Controller originals: `/mnt/data/slurm-results/wes-pilot/`; workers expose the
  same data as `/results/wes-pilot/`. Large immutable inputs/assets are pinned by
  their existing manifest hashes rather than copied into Git.

Branch: `feature/rust-wes-pilot-20261103`, based on
`2bee82a4d30ae298343e4045f50854e5feb36dc3`. Frozen deployments and helper binary
hashes are in each campaign manifest. The original Python-adaptive worktree
and its branch are not switched by this work.
