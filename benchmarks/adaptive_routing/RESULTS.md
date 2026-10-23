# Adaptive cache routing: engineering validation

2026-10-23 · `feature/adaptive-cache-routing`, based on `feature/fused-routing`
at `68ce69c18`. This is a new engineering experiment, not a replacement of the
publication campaign. Existing figures and manuscript measurements are untouched.

## Outcome

Adaptive BCF/CSI routing improves the cached path, but does **not** yet establish
a compelling WES advantage over a fast annotator. Keep it opt-in; legacy remains
the default. All 105 timed CLI conditions passed the separate complete-record
and relevant-header comparator, with no annotation-field exceptions.

| Real SGDP input | Legacy cached | Adaptive cached | Direct fastVEP |
| --- | ---: | ---: | ---: |
| WGS, LP6005519-DNA_A06 | 222.39 s | 147.06 s | Not rerun |
| WGS, LP6005592-DNA_C05 | 319.36 s | 240.49 s | Not rerun |
| Panel, 281 records | 26.24 s | 13.10 s | 11.74 s |
| WES, 82,714 records | 57.76 s | 52.31 s | 52.39 s |

WGS improvements are **1.51× and 1.33× over legacy cached**, not over direct
fastVEP. The real panel/WES are independent assay subsets of LP6005519-DNA_A06,
using the campaign's frozen ACMG/Twist regions. Each real condition was run once
with 16 allocated CPUs; subsecond differences are not evidence of an advantage.
Panel remains slower than direct fastVEP; WES is effectively at parity.

### Why high WES hit rates are not enough

The WES cache hit rate was 93.43% (5,432 misses). The single global annotation
recipe on those misses still took **37.22 s**, leaving **15.09 s** for the rest of
the timed CLI. Direct annotation took 52.39 s. Thus skipping 93% of records did
not skip 93% of annotator runtime. Startup/resource loading, miss composition and
other fixed work must be measured separately before assigning a cause; this
experiment does not isolate them. Approximately 10.46 s of the fused engine's WES
critical path was routing (see phase timings), so lookup is also unfinished work.

These runs tested **fastVEP, not vep-rs**. Beating direct vep-rs on WES remains a
target, not a demonstrated result or a guarantee. Next discriminating experiment:
matched direct-all/direct-misses/minimal-input runs for the same recipe and
resources, alongside cache-side phase profiling. If resource loading dominates,
further routing optimization alone will not produce large end-to-end gains;
persistent annotator integration or resource reuse becomes a separate direction.

## What changed

- Preserve explicit cache selection even with one worker; never silently scan
  the entire cache because a sample has only one contig.
- Generate an indexed normalized spool and 16-KiB query windows in one pass.
  Union conservative CSI chunk estimates; choose sparse windows below half the
  contig footprint, otherwise contiguous access.
- Schedule coordinate-disjoint tasks by estimated cache work plus sample work;
  split large contigs without splitting a POS group. Dynamic execution, ordered
  gather, one global annotation recipe invocation. Within dense tiles, work
  estimates remain approximate rather than perfectly cache-balanced.
- Pass decoder threads through; allocate foreground/codec work within the CPU
  budget and inherited affinity. Batch very small inputs to avoid process and
  filesystem startup overhead. Record strategies, allocation and phase timings.
- Require native helper 0.2.0; preserve existing BCF/CSI caches, immutable recipes,
  INFO overlay semantics, FORMAT values, sorted output and explicit fallbacks.

## Controls and mechanism

Three repetitions on the same 8-CPU worker, medians of CLI-only wall time:

| Cache-derived control | Legacy | Original fused (`c878d86c3`) | Adaptive |
| --- | ---: | ---: | ---: |
| 250 records, 22 contigs | 32.05 s | 45.32 s | 2.37 s |
| 250 records, one contig | 31.26 s | 110.64 s | 2.26 s |
| 75,000 records, 22 contigs | 34.49 s | 49.49 s | 10.21 s |

These controls are clustered, cache-derived, perfect-hit engineering inputs;
**they are not real WES samples or biological hit-rate estimates**. Additional
1-CPU/16-CPU cells exercise correctness and scheduling, not linear scaling.

An **untimed** strace control recorded cache BCF payload read calls falling from
1,026,089 (3,953,638,522 bytes) to 407 (1,821,047 bytes). These are userspace
file-read bytes, not physical disk traffic; CSI/index reads are excluded. Do not
use strace wall time as a performance result.

A final local dense synthetic control (one million records, mock annotator,
14 allowed CPUs, three alternating-order pairs) had engine-only medians
4.42 s original fused versus 3.89 s adaptive. The first pair passed the complete
comparator. This is not CLI time, real annotation or a real cohort speedup.

## Experimental mmap store: rejected for production integration

The explicit read-only prototype preserves raw BCF records and locus candidate
semantics. It is **not** an `annotate` backend and is never automatically built.

- Source cache: 18,869,857 records, 1,926,644,274 bytes.
- Derived store: 40,126,836,973 bytes (~20.8× source), predominantly raw payloads.
- Conversion: 503.96–539.73 s, reported separately and never reset on reuse.
- Median-hit WGS: 260.18 s store versus 147.06 s adaptive BCF.
- Lower-hit WGS: 334.92 s store versus 240.49 s adaptive BCF.

It failed both integration gates: ≥2× cache-handling and ≥1.5× end-to-end benefit
over **adaptive BCF**, not the old unrestricted scan. There is no amortization
payback because each run is slower. Per-run normalization, candidate BCF
materialization and indexing are included in timing. Raw payloads exceed a
worker's RAM and add substantial I/O; this rejects this prototype, not every
possible indexed format. Remaining candidates include avoiding hit-payload
spool materialization and reusing index/header handles, subject to measurement.

## Evidence, timing and limitations

- [Raw archive](../../benchmarks/results/raw/archive/vcfcache-adaptive-routing-20261023/README.md):
  frozen measured sources A/B/C/F, manifests, executable hashes, resources,
  complete comparisons, statistics and logs. Early native labels were 0.1.0;
  source/hash snapshots distinguish them. Final real-assay validation uses 0.2.0.
- [Final summary](final-summary-20261023/conditions.tsv) includes 93 control and
  12 real CLI conditions, grouped by experiment/sample/CPU/arm; distinct source
  snapshots and hardware are not pooled. The paired dense timer is separate.
- Separate Slurm steps per condition, CPU binding, light statistics; semantic
  audits outside timing. Use `resources.json:wall_seconds`, not the early outer
  timer that included Slurm startup. Largest-child RSS is not aggregate memory;
  use per-step cgroup peaks where available. CSI estimates are not observed I/O.
- Jobs 1257–1259 controls, 1261–1262 real WGS/store, 1265 trace, 1268 real assays
  finished. First assay setup 1266 failed before timed cells; 1267 was cancelled
  and replaced by 1268. An unused allocation on worker 5 had a credential error;
  no results from it are used. Queue was empty at final inspection.
- WGS and control results precede the final small-input/thread-allocation policy;
  frozen snapshots make this explicit. Final dense and real-assay validations
  exercise that policy. Do not claim all pilot cells measure the final binary.
- Large BCFs, derived stores and raw traces remain on the cluster, not in Git.
  This branch does not change published caches, release defaults or figures.

Recompute the CLI summary into a **new** output directory:

```sh
python benchmarks/adaptive_routing/summarize.py \
  --archive benchmarks/results/raw/archive/vcfcache-adaptive-routing-20261023 \
  --output /path/to/new-summary
```

Verification: local and Docker suites each **480 passed, 18 skipped**; native
unit tests **2 passed**. Network-gated tests remain skipped. Formatting, Ruff and
native Clippy checks are separate from semantic/performance evidence.
