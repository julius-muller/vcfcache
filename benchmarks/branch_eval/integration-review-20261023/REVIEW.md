# Adaptive integration review

2026-10-23. Branch: `feature/adaptive-integration`, based on
`feature/native-bcf-routing` at `224ca0fe8` (engine implementation `63757a07`).
This is a release candidate branch, **not a merge or a default change**.

## Decision

- **Adaptive is the production candidate.** It removes the old small-input
  traversal pathology, retains bcftools overlay semantics, and has broad wins
  with less bespoke typed-record machinery than resident native.
- **Native remains opt-in.** The paired fastVEP WES control below confirms a
  useful resident gain. Its additional Rust INFO/BCF handling and memory/IPC
  contracts do not yet justify making it the universal path.
- **Legacy remains the default and compatibility fallback.** The original
  full-cache-scanning fused generation is historical, not a release candidate.
  Current `fused` aliases adaptive; deleting it would break working callers.
- **No additional storage format or annotator integration.** The unsuccessful
  derived-store prototype is feature-gated out of default helper builds. Its
  source and archived evidence remain reproducible with `experimental-store`.

## Implemented maintenance work

Shared helper discovery/version validation, header validation, process handling,
affinity restoration and output staging now live in `engine_runtime.py`.
Public `adaptive` and opt-in availability-based `auto` engine names are added.
Auto uses adaptive with a supported optional helper, otherwise legacy if the
platform/helper is absent; explicit misconfiguration, stale helpers and runtime
errors fail. Native is never auto-selected. Stats preserve the historical
`engine` field while exposing canonical `effective_engine`, original request
and fallback/selection reasons. Existing cache layouts and recipes are unchanged.

The [backend-aware table](cells_backend_audit.tsv) preserves all 142 rows and
explicit supersession. Six native requests actually executed adaptive,
including VEP WES and every WGS/heavy request. Original tables/figures were not
overwritten. See [BACKEND_AUDIT.md](BACKEND_AUDIT.md).

## Focused confirmations, not a replacement campaign

### fastVEP WES: job 1327, sl-w2, 16 CPUs

Frozen deployment `63757a07`, same recipe/cache/input as job 1300. Two
alternating pairs; existing direct output reused for untimed complete-record
comparison. Every audit passed, with resident native confirmed in stats.

| Order | Backend | Wall seconds |
| ---: | --- | ---: |
| 1 | adaptive | 46.409 |
| 2 | native | 35.371 |
| 3 | native | 34.466 |
| 4 | adaptive | 43.040 |

Pair-wise native reductions: 23.8% and 19.9%; mean reduction ~21.9%.
These four cells support a resident benefit for this one case, not a universal
performance guarantee or a new cohort confidence interval.

### VEP WGS: job 1328, sl-w2, 16 CPUs

Reversed original order: `legacy@head` first, `main` second. Both reuse the
archived direct output for full-record audits. No expensive uncached run is
repeated. This job is **running** at the report snapshot; conclusions remain
pending. Results live at `/results/branch-integration/confirm-vep-1328` on
workers (controller: `/mnt/data/slurm-results/branch-integration/confirm-vep-1328`).
Do not treat submission as a completed comparison.

The Python ASTs of `_run_annotate` and `_run_annotate_nocache` are identical
between `136bafa9` and `63757a07`. SHA-256 of `ast.dump(...,
include_attributes=False)` is respectively
`7fae3f772e14613d8effece32ecf296977b13413de2efff7d8d11b61acdbd7df` and
`8e0da6878667bdb79b1083c71998ce698c2212b532e14660ca2ba750b7566194`.
This checks those two workflow methods, not the whole deployment or timing
environment; it is not evidence that every end-to-end cost must be identical.

## Clinical discrepancy: now independently reproducible

All 22 WES discrepancies are present in the original cache's CSQ payload.
Fresh annotation of the exact blueprint records agrees with direct sample
annotation. Removing GT/INFO/ID, repeating the run and using one thread do not
recreate the cache payload. Blueprint neighborhoods of 1 kb still agree with
direct; 100 kb (3,110 input records), and full chr1 (1,435,577 records), recreate
all 22 cached differences. The binary and ClinVar data/index hashes match cache
build provenance. Transcript-tree identity was recorded at build time but has
not been independently rehashed during these focused probes.

Delta debugging gives this sorted, three-record input:

| CHROM | POS | REF | ALT | Role |
| --- | ---: | --- | --- | --- |
| chr1 | 235699670 | G | T | preceding context |
| chr1 | 235700300 | T | C | preceding context |
| chr1 | 235782202 | A | G | target |

With all three records, the target's `dbNSFP_CV_CLNSIG`, `dbNSFP_CV_REVSTAT`
and `dbNSFP_CV_CLNDN` values are empty. Target alone, or with either one
preceding record, yields `Benign`, `criteria_provided&_single_submitter`,
`not_provided` for all four target consequences. This reproduces at `--fork 1`
and `--fork 16`, both with the original recipe and **dbNSFP alone** (gnomADc
removed). No VCFcache engine is involved.

This localizes the demonstrated omission to the tested vep-rs clinical recipe's
dbNSFP lookup/annotation path and its input context. It does **not** establish
the precise source-code mechanism, a defect in every vep-rs release, or that
every one of the 40 WGS discrepancies has that same minimal trigger. The
[upstream v0.4.0 dbNSFP implementation](https://github.com/natera-open-source/vep-rs/blob/v0.4.0/crates/vep-builtins/src/plugins/dbnsfp.rs)
uses a batched annotation store; source review alone does not prove binary
identity or isolate the failing internal component.

Evidence: [minimal input](minimal-clinical-1331/minimal.vcf),
[reduction](minimal-clinical-1331/result.json),
[thread/plugin controls](clinical-confirm-1332/results.json), and the
`clinical-1326`/`clinical-1329` command/output snapshots. The reproducer retains
original absolute resource paths; adapting it to another host requires those
same pinned resources. No upstream issue, message or PR was sent.

**Publication status:** clinical speedups remain provisional timing evidence,
not output-equivalent validated acceleration. No fields are ignored, cached
payloads repaired in place, or failed comparisons relabelled as passes. An
upstream fix/workaround needs a new immutable recipe/cache plus full WES/WGS
audits before making that claim. Core vep-rs and fastVEP validity are separate.

## Release gates and verification

- Host suite: **547 passed, 21 skipped**; container Python 3.13 suite: **547
  passed, 21 skipped**, using current source and rebuilt Linux helper 0.4.0.
- Rust unit tests, default-feature clippy, formatting and experimental-store
  compilation passed; Python ruff/black passed.
- Tests cover explicit aliases, automatic optional-helper/platform fallback,
  invalid/stale helpers, exact-record parity, native fallback, affinity cleanup
  and staging cancellation. Existing stdout/BCF/header/allele tests remain.
- Default helper excludes store commands. Normal Python/lean-image runtime
  remains Rust-free. Default change still requires release packaging/platform
  support and representative integration smoke checks; current default remains
  legacy.
- A fresh complete Docker image build hit container DNS failures during pip
  installation. The Rust Docker stage built successfully; the full container
  suite used the existing dependency image with current source and that new
  helper mounted read-only. This is not a claimed fresh image-build success.

## Separate next experiment: cohort miss batching

No more engine rewrite is justified at the heavy end. Before a public batching
API, run a small isolated control with a **variant-local** recipe:

1. Gather per-sample misses, deduplicate full allele keys, annotate their union
   once, and overlay annotations back while preserving each sample's original
   GT/FORMAT and record metadata. Do not union sample-dependent annotations.
2. Compare against both per-sample direct annotation and **batched uncached
   union annotation**. Equalize CPU quotas, startup treatment and recipe
   resources; separate deduplication savings from invocation amortization.
3. Time union construction, deduplication, demultiplex/overlay, compression and
   indexing. Report setup/build costs separately; strict per-sample full-record
   audits remain outside timed cells. Record effective backends and resource
   hashes. Reusing original direct outputs is legitimate only for unchanged
   recipes/resources and with documented context-invariance.
4. Start with a few public samples and one recipe, not another broad campaign.
   Do not use the failing clinical dbNSFP recipe until its context sensitivity
   is resolved. Only add a product API if throughput gains beat the fair direct
   batching baseline and maintain exact sample output.
