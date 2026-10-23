# Resident native BCF engine: results and limitations

2026-10-23 · branch `feature/native-bcf-routing`, based on adaptive routing
commit `8830bf7b3`. Engineering pilot only; existing publication figures and
manuscript measurements have not been replaced.

## Outcome

The resident engine improves cache handling without changing BCF/CSI, recipes,
or annotator interfaces. It gives a validated fastVEP WES benefit, but the
goal of beating **core vep-rs** on WES was **not reached**. The clinical vep-rs
comparison exposed annotation differences and is not an acceptable speedup claim.
Keep this prototype opt-in; legacy remains the default.

Primary runs use frozen deployment B, one run per real sample/condition,
16 allocated CPUs, local cache BCFs and the same recipe/resources in each arm.
Complete-record/relevant-header audits are outside timing, with no field
exceptions. Sample names below abbreviate `LP6005519-DNA_A06` and
`LP6005592-DNA_C05`. WES inputs are assay subsets of the existing SGDP WGS
files; this pilot makes no new relatedness or database-membership claim.

| Annotator | WES sample | Records | Hit rate | Adaptive fused | Native | Direct |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| fastVEP clinical | A06 | 82,714 | 93.43% | 51.14 s | 41.09 s | 53.05 s |
| fastVEP clinical | C05 | 108,539 | 71.12% | 54.61 s | 46.26 s | 55.25 s |
| vep-rs core | A06 | 82,714 | 93.43% | 23.71 s | 12.94 s | 9.42 s |
| vep-rs core | C05 | 108,539 | 71.12% | 17.82 s | 15.35 s | 9.32 s |

All four sample/recipe comparisons passed. For fastVEP, native saves **11.96 s
and 8.99 s versus direct** (1.29× and 1.19×). Core vep-rs remains faster uncached.
Do not pool measurements from different nodes or interpret single-run fractions
of a second as an advantage. Annotation references retain existing shared paths;
resource-load timing can vary even though every arm uses identical resources.

Deployment C then changed the final merge to consume and release records while
writing, instead of accumulating all recipe INFO payloads. Same-worker checks
against saved direct outputs passed: 42.91 s fastVEP and 13.46 s vep-rs core on
A06. These are revalidation cells, not independent biological replicates.
The final repository additionally fixes returned-record statistics for recipes
with missing required tags and conservatively replays case-insensitive matches.
Frozen sources and executable hashes distinguish these development snapshots.

## What the implementation removes

- One resident Rust process retains original WES records and cached INFO.
  bcftools splitting streams through a FIFO; no indexed normalized sample spool.
- Native HTSlib BCF/CSI readers query occupied 16-KiB windows. Unrelated cache
  candidates are not copied/unpacked into full Rust records. A shared codec pool
  and lookup workers stay within the configured CPU affinity/budget.
- One indexed miss BCF is annotated by the unchanged external recipe once.
  Exact INFO merging preserves original FORMAT/GT, identifiers and other fields.
- One compressed final output and CSI replace hit spools, shard assembly and
  compressed-shard gather. Final recipe merging holds only one output locus at
  a time while releasing originals.
- Buffer-budget/unsupported-input fallback happens before annotation. A rare
  compatibility merge reuses completed recipe output rather than rerunning it.
  Parsing or subprocess failures never silently become successful fallbacks.

This is not a new cache format or in-process annotator. There is no vep-rs API
integration, new service, custom key/value store or automatic mutable cache.
Remaining overhead is measurable; linear scaling is not assumed.

## Assay projection: useful for lookup, not a guaranteed end-to-end win

Only after full-cache measurements, the same annotation BCF was projected onto
the frozen Twist WES BED. The parent cache and recipe stayed unchanged. Projection
is ordinary BCF/CSI, records parent/BED identity and construction cost, and treats
out-of-projection positions as misses rather than dropping sample records.

| Annotator | Projected BCF | Construction | A06 native / direct | C05 native / direct |
| --- | ---: | ---: | ---: | ---: |
| fastVEP clinical | 87.12 MB | 70.07 s | 39.23 / 50.41 s | 50.01 / 53.45 s |
| vep-rs core | 39.98 MB | 54.10 s | 11.06 / 9.25 s | 13.08 / 9.00 s |

Sizes exclude CSI. Parent BCFs are 1,926.64 MB and 840.32 MB respectively.
Native lookup fell from **3.36 to 1.65 s** for fastVEP A06 and **2.68 to 1.50 s**
for core A06. Both projected samples have 75 additional misses because the input
region extraction includes overlapping records whose POS lies outside the BED;
the projection uses POS selection. All projected full-record audits passed.

The lower-hit fastVEP sample's total time increased despite faster lookup:
its recipe stage took 40.49 s versus 35.59 s in the full-cache condition.
Thus the small lookup benefit is not evidence of a reliable total improvement.
Do not recommend automatic projection or quote an amortized payoff from these
single-run timings. It is an experimental preparation utility, not a public alias.

## Clinical vep-rs: failed semantic gate, not a validated 3.6× speedup

On A06, fused and native agreed, but fresh direct annotation differed at
**22/82,714 records**, in ClinVar-derived CSQ values. Native took 18.85 s and
direct 67.50 s, but quoting their ratio as a package speedup would violate the
output-equivalence requirement. The runner stopped before a second sample or
clinical projection; no fields were ignored and no cache payload was rewritten.

Two discrepant loci were rerun in isolation. Those results matched full direct
annotation and still disagreed with cache-derived annotations. The evidence
isolates this from the native/fused routing difference; it does **not** establish
the underlying vep-rs/plugin cause. Input context, source-cache generation and
recipe cacheability need a separate investigation before using this clinical
recipe for publication or advertising. Resource/recipe pinning alone does not
prove that every annotation is variant-key independent.

New vep-rs caches were locally built from the bundled AF>=1% blueprint, **not**
downloadable bundled vep-rs caches. True one-time build durations are retained:
130.43 s core, 778.47 s clinical. The full transcript conversion completed
3,102/3,102 files, 552,279 transcripts without warnings. A redundant params-only
version-command key was removed in freshly staged recipe metadata; original
recipes, annotation commands and true build costs were preserved in provenance.

## Startup calibration and the performance target

Direct full-CLI calibration, same recipe and worker; the minimal inputs retain
one miss from each occupied contig. Preparation and audits are untimed:

| Recipe | 5,432 misses | 22 variants / same contigs |
| --- | ---: | ---: |
| fastVEP clinical | 39.00 s | 12.11 s |
| vep-rs core | 7.37 s | 6.61 s |
| vep-rs clinical | 13.98 s | 7.58 s |

These are input-size/composition diagnostics, **not pure startup measurements**.
The fastVEP minimal-input result does not support treating the two-point fitted
36-second intercept as universal per-invocation startup. Miss composition and
nonlinear annotation work must be accounted for. The chr22 slopes from job 1251
likewise are not measured whole-WES latencies.

The ≤2 s cache-side target was not reached. For fastVEP A06, full native CLI
minus its recipe stage was about 8.14 s; the projected case was about 6.20 s.
Reducing lookup alone therefore cannot deliver the complete remaining goal.
The next decision should use these phase timings, not assume a different file
format or annotator integration is needed. This iteration stops without expanding
into either, as requested.

## Controls, verification and provenance

The short repeated control uses genuine WES locus distribution but is deliberately
restricted to reusable variants and stripped of annotations before lookup. Its
perfect hit rate is an engineering control, **not a biological cohort claim**.
Three paired runs alternate arm order; their medians and all raw timings remain
separate from the single-run real samples: **12.79 s fused versus 6.03 s native**
(2.12× cached-path improvement; 77,282 records). All six full-record audits passed.
The final WGS test exercises bounded memory fallback, not a claim that the resident
backend accelerates WGS: **4,858,092 records** matched completely, including all
INFO/FORMAT values. Both requests used fused routing after the native request's
indexed-count budget check; their CLI times were 112.05 and 112.64 s.

All submitted work has finished; the queue was empty at final inspection.
The [complete summary](summary-20261023-complete/conditions.tsv) retains 44 timed
cells, including one successful baseline from an earlier failed diagnostic pilot.
The only failed cross-condition semantic gate is the explicitly flagged clinical
vep-rs comparison above. Baseline-versus-itself audits validate indexed readability,
not independent annotation correctness; the cross-arm audits provide that check.

Local and Docker suites: **525 passed, 18 skipped**. The native differential
module contains 45 cases covering complete records/headers, flags and missing
INFO, arbitrary INFO dictionary IDs, window-crossing deletions, shared positions,
missing contigs, dropped annotations, recipe ordering, compatibility reuse,
stdout, memory fallback and affinity restoration. Native unit tests, Clippy,
Ruff, formatting and scoped typing checks are separate release gates.

An early pilot's BCF read sequentially but failed indexed access. The fix flushes
the BGZF header before initializing on-write CSI and checks retrieval before
publication. A large-header regression test covers it. Failed setup attempts and
the failed clinical audit are retained, not hidden or relabelled successful.

[Raw evidence and frozen sources](../results/raw/archive/vcfcache-native-bcf-20261023/README.md)
include manifests, comparisons, resources, original build costs and setup logs.
Large BCFs/transcript resources stay on the cluster. No existing figures,
publication assets, cache aliases or release defaults were changed.
