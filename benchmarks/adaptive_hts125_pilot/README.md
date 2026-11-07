# Python-adaptive / HTSlib 1.25: bounded follow-up pilot

This branch is based on `feature/python-adaptive` at `8ec84894f`, not the newer
input-preserving main implementation. Nothing is promoted to main by this pilot.
No Rust/native engine, new cache format or full-cohort campaign is introduced.

## Questions separated by controls

| Condition | Source | bcftools/HTSlib | Input pass | Lookup / shards |
|---|---|---|---|---|
| baseline-121 | Frozen original Python-adaptive | 1.21 | Original | Original auto |
| baseline-125 | Same original source | 1.25 | Original | Original auto |
| prepared-auto125 | Pilot branch | 1.25 | None | Original auto/budget |
| prepared-contigs2 | Pilot branch | 1.25 | None | Whole contigs / 2 requested |
| prepared-windows2 | Pilot branch | 1.25 | None | Merged 16-kb windows / 2 |
| prepared-windows4 | Pilot branch | 1.25 | None | Merged 16-kb windows / 4 |
| direct-125 | Pilot branch, uncached | 1.25 | None | One annotator invocation |

The first pair isolates the external toolchain. Original-auto versus prepared-
auto tests removing the redundant input pass (exact matching is also enabled;
the prepared domain makes this representation-preserving). Contigs2 versus
windows2 isolates access geometry; windows2 versus windows4 compares concurrency.
Prepared-auto versus contigs2 also changes shard count/codec budgeting and must
not be attributed to only one mechanism.

## Verified experimental domain, not a generic product fast path

`VCFCACHE_PILOT_PREPARED=1` is a private trusted-input assertion. Before timing,
every source is checked to be biallelic, free of ALT=*, unique by variant key,
and free of existing recipe annotations. Unsupported inputs fail instead of
being normalized. Both candidate direct and cached runs reuse that same indexed
BCF. Original input/cache files are never changed. Assay extraction/encoding,
asset hashes and input checks are untimed; region construction and normal
lookup/output work remain timed. Full record/header comparisons occur afterward
without ignored fields. No new full-file recipe fingerprint scan is added to
the timed path. Existing cheap completeness/index/pipeline checks remain.

`VCFCACHE_PILOT_LOOKUP` chooses auto/contigs/windows only for this experiment.
`VCFCACHE_PILOT_MERGE_WINDOWS=1` coalesces adjacent occupied 16-kb blocks.
`VCFCACHE_PILOT_THREAD_BUDGET=1` reserves one main process per shard and splits
remaining codec workers across the requested 2/4 shards. Effective shard count
and codec budget are recorded. All timed steps have enforced eight-CPU affinity.
Wrappers thread supported commands inside the frozen recipes, but never
override explicit per-shard thread options. Codec features stay matched to the
existing zlib-based 1.21 build; no compression-library experiment is added.

## Campaign and later review

Five exclusive blocks occupy Intel Broadwell workers **sl-w1/2/3/4/6**;
**sl-w5 (AMD) is excluded**, including from collection. Same KPGP-00139 inputs
and frozen recipes/caches as the preceding pilot: vep-rs core WES/WGS, fastVEP
clinical WES/WGS, and vep-rs clinical WES (diagnostic). One fresh direct per
block, one timing per condition, rotated condition order, no fastest-run
selection and no extra technical repeats. OS/file caches are not flushed.
Shared NFS activity can still add variability; same-worker pairs and phase data
are engineering evidence, not a biological-cohort uncertainty estimate.

Every cell is checkpointed with commands, source/binary/input/cache hashes,
CPU/storage metadata, counts, resources and comparisons. Five seven-condition
blocks mean **35 planned timings** and 30 cached-versus-direct comparisons.
Known clinical annotation mismatches are diagnostic, never validated gains.

After-any Slurm collection writes TSVs, a Markdown report and a frozen evidence
archive. A bounded controller postprocessor renders R/ggplot2 PDF/SVG/600-dpi
PNG figures in a new folder, including separate panels. Neither process
submits further annotation work. Failed/time-limited scratch is retained for
review days later. The safety wall limit is two hours; expected first WES
results are roughly 10 minutes after launch, with the full pilot/audits typically
30–60 minutes. These are estimates, not completion claims.
