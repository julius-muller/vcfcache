# Fixed-worker decoder-pool screen — 3 November 2026

## Outcome

**A modest lookup optimization, not a reliable end-to-end improvement.**
Keeping four lookup workers and disabling the extra decoder pool speeds cache
lookup in five of six paired WES cases. Miss-recipe time dominates the full CLI,
however, and total-runtime gains are inconsistent. Keep the new override
experimental; do **not** change the native default on this screen.

| Recipe | Inputs | Median paired lookup gain | Lookup improved | Median paired full-CLI gain | Full CLI improved |
| --- | ---: | ---: | ---: | ---: | ---: |
| fastVEP clinical | 3 | 1.52× | 3/3 | 1.02× | 2/3 |
| vep-rs core | 3 | 1.40× | 2/3 | 0.94× | 1/3 |

Ratios are four-decoder-thread time divided by zero-decoder-thread time, so
greater than one favours disabling the pool. These are medians of paired
ratios, not ratios of marginal medians. Median measured lookup savings are
1.62 s for fastVEP and 0.96 s for core. Median full-CLI savings are 0.90 s and
−1.05 s, respectively.

## All paired observations

Each cell is a single timed invocation; no best-run selection or technical
repeat averaging is used. Both conditions are retimed here, unlike comparisons
to the older pilot's cached controls.

| Recipe | Input | Lookup, 4 decoders → 0 (s) | Full CLI, 4 decoders → 0 (s) |
| --- | --- | ---: | ---: |
| fastVEP clinical | KPGP-00139 | 4.74 → 3.13 | 46.05 → 46.46 |
| fastVEP clinical | LP6005592-DNA_C05 | 4.93 → 3.07 | 63.41 → 54.81 |
| fastVEP clinical | KPGP-00231 | 4.81 → 3.48 | 43.20 → 42.30 |
| vep-rs core | KPGP-00139 | 4.15 → 5.18 | 15.61 → 21.99 |
| vep-rs core | LP6005592-DNA_C05 | 3.34 → 2.38 | 17.94 → 17.24 |
| vep-rs core | KPGP-00231 | 4.24 → 2.60 | 16.71 → 17.76 |

The largest fastVEP total-time gain is **not** 8.61 s of decoder savings:
only 1.86 s is measured in lookup, while the miss-recipe invocation changes
from 52.76 to 45.90 s. Conversely, core KPGP-00231 saves 1.64 s in lookup but
its recipe invocation increases from 8.62 to 11.18 s, leaving the full CLI
slower. Single-run recipe differences must not be attributed wholesale to the
decoder switch. This screen is too small for a universal dispatch rule.

## What was isolated, and what was not

- Four lookup workers in both arms; the same eight-CPU affinity allocation,
  unchanged annotator/output threading, standard BCF caches, recipes, input
  preparation, 16-kb query windows, and 2-MiB block-cache setting.
- Explicit zero versus four shared cache-codec threads. No changes to lookup
  worker count, coalescing, whole-cache residency, or annotation settings.
- Six exclusive Slurm jobs, one per sample/recipe pair, running across all six
  workers. Each condition has a separate fresh Slurm step and process. Order is
  counterbalanced across the six pairs: three start with each setting. Within
  either recipe the three-input allocation necessarily splits 2:1, not 1:1.
- Filesystem/page caches are not flushed. This is application-start latency,
  not a controlled disk-cold experiment. No technical replicates were added.
- Original measured direct baselines are reused only for context in
  `conditions.tsv`; the decoder decision uses the newly paired cached arms,
  not a comparison against an old cached timing.
- Disabling the pool jointly removes asynchronous decode, prefetch and
  seek synchronization, while allowing the ordinary decoded-block cache to
  populate. This does **not** isolate cache hit rate from threading overhead.

Source inspection supports the mechanism but not a universal speedup: in the
bundled HTSlib, `bgzf_mt_seek` resets the decode queue and discards queued work.
The multithreaded reader probes `load_block_from_cache`, but only the
single-threaded read path calls `cache_block` to insert decoded blocks. New
pooled readers therefore cannot build the same ordinary block-cache reuse.
The earlier shorthand that the threaded reader simply never probes the cache
was incomplete. A larger block cache was not tested, and no HTSlib patch was
needed for this experiment.

## Implementation and validation

Helper **0.8.4** adds an optional `codec_threads` job field. If absent, the
existing `CPU budget − effective lookup workers` rule is unchanged. Opt in with
`VCFCACHE_EXPERIMENT_CODEC_THREADS=0`; the controller requires helper >=0.8.4.
Both controller and helper reject allocations exceeding the available CPU
budget. `lookup_codec_threads` records the effective setting, and the pilot
audits that it actually matches the requested condition.

All **12 timed outputs pass strict complete-record/relevant-header comparison**
against their retained direct reference, with no ignored fields. Within each
pair, input/miss/output counts, candidate-record counts, query counts and
semantic SHA-256 digests agree. No native fallback is accepted. Nine new tests
cover pool settings, zero/partial/full hits, and invalid budgets; the full suite
passes **711 tests, 21 skipped**. Ruff, Black and Rust formatting checks pass.

## Recommendation

Retain the small opt-in control for investigation; do not disable the decoder
pool universally. It can save roughly a second or two in these sparse lookups,
but this does not resolve core-WES break-even or demonstrate a new 2× full-run
gain. A future region-driven policy would need a larger paired screen and a
dense-scan control before choosing decode mode automatically. More lookup
optimization alone cannot remove the measured miss-recipe fraction.

This experiment is finished. No further jobs, automatic follow-up, publication
rerun, or default change is queued.

## Provenance and files

- Source: `f72a8f9916e6b42d70328c1d883a5dcbe0396166` on
  `feature/rust-wes-pilot-20261103`.
- Helper SHA-256:
  `53a42e29f25590ac5df6ef2e149fc2322a7f5ea30a7c3cadacb6aa49584bcef8`.
- Jobs **2743–2748**, all completed successfully.
- Controller: `/mnt/data/slurm-results/wes-pilot/20261103-decoder-screen`;
  worker alias: `/results/wes-pilot/20261103-decoder-screen`.
- `conditions.tsv`, `pairs.tsv`, `summary.tsv`, `attempts.tsv`: all observations,
  paired effects and acceptance states, including regressions.
- `terminal-evidence.tar.gz`, `EVIDENCE_INVENTORY.json`, `SNAPSHOT.sha256`:
  manifests/READY, submitted commands, resources, statistics, strict comparisons
  and hashes. Extract `source/` here and run
  `python -m benchmarks.wes_pilot.decoder_collect benchmarks/wes_pilot/decoder-results-20261103`
  to regenerate the tables. Direct-reference receipts are retained in the
  earlier pilot's immutable terminal-evidence archive.

Earlier pilot results and printed figures are unchanged. No figure was added
for this six-pair diagnostic screen.
