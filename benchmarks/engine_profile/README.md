# Bounded native/adaptive overhead investigation

Engineering profiling only. This branch does not change the default engine,
annotation semantics, cache format or manuscript. Baseline: `296855d72`.

## Frozen design

86 cells (cap 96), including 24 once-only reference cells and 24 primary cells.
Two independent SGDP genomes (`LP6005519-DNA_A06`, `LP6005592-DNA_C05`), their
BED-derived WES projections, and fastVEP clinical, vep-rs core and VEP everything.
The context-sensitive vep-rs clinical recipe is deliberately excluded.

Primary conditions run once, with native/adaptive order reversed for the second
sample. Focused diagnostics use the first sample: 1/8/16 enforced CPUs, 2/4/8
workers, native 512/2048/8192 MiB buffer budgets, BCF versus VCF.gz, and skipping
splitting only after verifying biallelic input. Three alternating instrumentation
on/off pairs quantify disturbance. Two optional `perf` cells explicitly skip if
sampling permission is unavailable. First-use-after-copy and repeated-access do
**not** represent controlled cold/warm page-cache measurements. No caches are
dropped. Staging and fingerprinting warm files.

Fresh direct and uninstrumented cached references are generated once for these
new frozen fixtures; historical publication results are neither overwritten nor
silently reused without an input/output identity binding. Every measurement is
compared to both references outside timing. Every INFO/FORMAT field and relevant
header is checked. A known difference versus direct annotation stays a failed
direct audit even if the profiling-equivalence audit passes. Nothing is ignored.

## Deployment and preparation

Deploy two immutable source snapshots and compatible compiled helpers: baseline
and this profiling commit. The request JSON has `builds.baseline` and
`builds.profile`, each with `root`, `commit`, `files` (relative-path SHA256 map),
and helper version/build metadata. Both use the same frozen Python environment.
Run on a worker, as appuser:

```sh
PYTHONPATH=PROFILE_DEPLOYMENT PYTHON -m benchmarks.engine_profile.run prepare \
  /results/engine-profile/UNIQUE_CAMPAIGN --request REQUEST.json
```

Preparation claims an empty campaign, requires 300 GiB shared free space, copies
caches and sample fixtures, freezes binaries, recipes, FASTA/indexes, transcript
and supplementary databases, and writes `manifest.json` plus its SHA256 in
`READY`. The original inputs and caches are not edited. An obsolete vep-rs recipe
metadata key, if present, is removed only from the new copy; original bytes are
retained as `annotation.original.yaml`. Each group verifies frozen resources and
code again before measurement and requires 60 GiB local scratch headroom.

## Submit and leave

On the controller, shared storage is `/mnt/data/slurm-results` instead of
`/results`; use the controller path for CAMPAIGN and PYTHONPATH:

```sh
PYTHONPATH=PROFILE_DEPLOYMENT PYTHON -m benchmarks.engine_profile.run submit CAMPAIGN
PYTHONPATH=PROFILE_DEPLOYMENT PYTHON -m benchmarks.engine_profile.run status CAMPAIGN
```

Six pinned-node groups request 16 CPUs, 28000 MiB, exclusivity and six hours each.
Two groups share sl-w2 and run sequentially. sl-w5 is intentionally excluded.
Each cell gets a **new srun step/cgroup**, enforced CPU affinity, 1-second durable
resource samples, commands and logs. CPU/I/O deltas and independent memory peaks
are retained. An `afterany` collector is registered even after partial submission.
Accepted job IDs are written immediately to `submissions.json`. No SSH session,
terminal multiplexer or submitting Python process must remain alive.

Reference and failed-audit output BCFs/indices are retained on shared storage;
successful measurement outputs can be reconstructed from immutable inputs and
commands. Small logs and metadata are copied before owned per-cell scratch is
removed. Failed scratch is retained. Group-level exceptions are recorded.

## Review days later and explicit retries

```sh
PYTHONPATH=PROFILE_DEPLOYMENT PYTHON -m benchmarks.engine_profile.run collect CAMPAIGN
PYTHONPATH=PROFILE_DEPLOYMENT PYTHON -m benchmarks.engine_profile.run resume CAMPAIGN \
  --cells g00-c012 g00-c013
```

`status` is read-only. `collect` retains terminal scheduler evidence before its
retention expires and creates content-addressed reports; it also works offline
after copying the campaign. Missing/running work is never counted as successful.
Retries require explicit cell IDs, refuse successful cells and active prior jobs,
and write a new numbered attempt. Resume failed references and dependent skipped
cells together. The paired storage-access experiment requires both cells in the
same attempt; retrying a lone repeated-access cell is not a valid replication.

`latest-review.json` points to `REPORT.md`, `cells.tsv`, `phases.tsv`,
`instrumentation.tsv` and `summary.json`. JSONL traces carry a shared epoch clock
and monotonic durations, including abandoned native attempts before fallback.
Inclusive spans overlap: **do not sum phases or worker durations as wall time**.
Only interval unions are reported as covered time. Coarse spans cannot separate
compression/decoding CPU; use permitted stack sampling or retain that uncertainty.

Plots are generated exclusively with R/ggplot2:

```sh
Rscript benchmarks/engine_profile/plot.R REPORT_DIRECTORY
```

If R is absent on the collector, tables and plot instructions still survive.
Collection never overwrites prior manuscript figures. A candidate optimization
must plausibly save at least 5% end-to-end, or address meaningful fallback/resource
waste, with acceptable complexity and strict equivalence. Phase elimination is
an optimistic upper bound, not an implemented or promised speedup. Check the
three on/off pairs before interpreting small timing differences. This campaign
does not itself implement any optimization or select a production backend.
