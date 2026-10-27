# Performance decision: release candidate versus another optimization round

27 October 2026. Engineering assessment of the frozen 23 October profiling
campaign and the completed 25/27 October VEP WGS recovery. No engine change,
new benchmark submission, default switch or manuscript overwrite is part of this
consolidation.

## Recommendation

**The product story is already convincing for annotation-rich workflows and
fastVEP WGS. It is not yet a convincing universal acceleration story for minimal
vep-rs, especially WES.** Do one bounded optimization/validation sprint before
freezing the release candidate and rerunning the publication series. Do not
start an open-ended rewrite, and do not rerun the whole series while engine
selection and the rich vep-rs recipe are still undecided.

Highest priorities:

1. Obtain a strictly equivalent, genuinely annotation-rich vep-rs recipe. This
   has greater relevance to the desired headline than making an already cheap
   consequence-only recipe marginally faster. Existing clinical timings are
   promising but cannot yet substantiate the validated headline.
2. Reduce the measured routing/assembly overhead, starting with low-risk
   profiling-guided experiments. Target **2x rather than 1.4–1.5x core-vep-rs
   WGS**, with the lower-hit genome retained as a gate. This is a target, not a
   forecast.
3. Keep adaptive as the broadly applicable production candidate and native as
   a useful opt-in WES path. Retain the existing common orchestration rather
   than growing four independent implementations. Legacy is still the current
   default; changing it requires an explicit release decision.

## What was consolidated

- 86 original planned conditions; eight VEP WGS conditions recovered separately
  after two direct-reference OOM failures and a later group timeout.
- The final 86-condition view contains **76 strict-pass completions, eight VEP
  WES completions with direct-output differences, and two optional perf skips**.
  All 24 primary measurements executed. The original failed and timed-out
  attempts remain in the evidence; they are not silently erased.
- All eight recovery cells passed. Four VEP WGS primary requests actually ran
  adaptive: their native requests fell back before resident annotation.
- Both genomes are SGDP: A06 has 82,714 WES / 4,858,092 WGS records; C05 has
  108,539 / 6,169,350. Observed cache hit rates are **93.4% / 92.7%** for A06
  WES/WGS and **71.1% / 79.0%** for C05. These are observed inputs, not universal
  WES/WGS expectations or proof about population relatedness.
- VEP WGS recovery uses four annotation forks and a 50,000-record buffer;
  original VEP WES uses 16 / 500,000. Engine CPU allocation remains 16. The
  configurations are never pooled as one per-variant scaling experiment.

## User-facing acceleration already demonstrated

The table uses **uninstrumented adaptive cached references**, not the fastest
arm selected after looking at results. Each comparison uses its own sample,
recipe and direct baseline. Time savings exclude cache construction.

| Recipe / assay | A06 direct → cached | Speedup | C05 direct → cached | Speedup | Validity |
| --- | ---: | ---: | ---: | ---: | --- |
| fastVEP clinical WGS | 446.1 → 111.1 s | **4.01x** | 665.8 → 230.3 s | **2.89x** | Strict pass |
| vep-rs core WGS | 66.9 → 46.1 s | **1.45x** | 82.7 → 60.9 s | **1.36x** | Strict pass |
| VEP everything WGS, recovery | 4.07 h → 21.4 min | **11.41x** | 5.73 h → 65.4 min | **5.25x** | Strict pass |
| fastVEP clinical WES | 47.3 → 47.3 s | 1.00x | 62.3 → 54.2 s | 1.15x | Strict pass |
| vep-rs core WES | 10.85 → 17.00 s | 0.64x | 9.29 → 19.07 s | 0.49x | Strict pass; slower |
| VEP everything WES | 268.5 → 133.7 s | 2.01x | 321.8 → 178.7 s | 1.80x | Two differing records per sample |

VEP WGS saves **3.71 and 4.64 hours/sample**. FastVEP WGS saves approximately
**75% and 65%** of direct runtime, including a genome with 79% cache hits. This
is a strong, intuitive benefit, even with a fast annotator. Neither the core WES
loss nor the clinical equivalence caveat should be hidden behind these examples.

Resident native improves the fastVEP WES primary measurements to **1.18x and
1.19x** versus direct, and core-vep-rs WES to **0.81x and 0.52x**. Native's
incremental gain matters, but is sample-dependent: compare these with the
paired alternating fastVEP WES control in the earlier integration review,
which found a 20–24% improvement over adaptive on its one case. Do not advertise
a universal native percentage or call native fallback a native speedup.

The old branch matrix also contains substantial panel improvements, but this
profiling campaign did not run panels. Preserve those measurements as separate
historical evidence rather than inventing new panel confidence from WES/WGS.

Historical annotation-rich vep-rs timings reach **3.42x WES** and approximately
**3.7x WGS**, but the pinned clinical recipe has 22/40 ClinVar-derived CSQ
differences. A three-record reproduction works without VCFcache. That explains
why routing is not implicated; it does **not** waive output equivalence. A fixed
or different real recipe must pass full-record audits before those ratios can
become publication claims. This assessment does not claim whether a newer
upstream release has fixed it.

## Where there is still leverage

### 1. Core-vep-rs WGS: a real routing and reconstruction opportunity

The A06 adaptive primary takes **43.80 s**, of which **8.75 s** is the measured
miss-annotation stage. The remaining **35.04 s** is cache-side orchestration,
preprocessing, lookup, output construction and other work. Its route interval
is 12.10 s and assembly interval 10.16 s; preprocessing and final indexing each
take about 3 s. These are interval unions, not sums of parallel workers.

C05 takes **59.87 s**, with **16.02 s** of miss annotation and **43.85 s** of
remaining work. With unchanged miss annotation, reaching 2x overall requires
reducing those residuals by **29.5% and 42.2%**, respectively. Reaching 3x would
require **61.3% and 73.7%** reductions. The latter is a much more ambitious
architectural target, not a small tuning exercise.

Inspect compressed/uncompressed intermediate materialization, rereading miss
records, assembly overlay, and the final index pass. Preserve full-record,
allele-pairing and header semantics. Current code already uses uncompressed
intermediates and naive BGZF concatenation in places; do not propose them as
unimplemented discoveries. The [bcftools manual](https://www.htslib.org/doc/bcftools.html)
documents streaming support but also indexing requirements for multi-input
operations. Removing materialization is an experiment with correctness and
access-pattern constraints, not automatically a win.

### 2. Native WES: smaller, targeted lookup and startup improvements

A06 native looks through **2.97 million candidate cache records** for 82,714
input records, using 14,147 queries. Native source rounds requests to 16-kb
windows. This is roughly 36 candidates per input record, **not 36x redundant
physical I/O**: block granularity, shared decompression and seeks matter.
Lookup costs 2.61 s for core and 3.35 s for fastVEP. A tighter window/block
selection experiment is worth trying; simply issuing point queries could be
slower. Candidate counts alone cannot predict the gain.

Core WES native has 6.47 s of residual work on A06 and 10.77 s on C05. With the
observed annotation stage unchanged, parity requires reducing these residuals
by **39% and 80%**. Even deleting *all* residual work would yield only
approximately **1.57x and 1.30x** direct speedups. Thus a stunning consequence-only
WES headline is not supported by this implementation's measured miss stage.
This is a conditional bound, not a theorem about all future designs.

C05 also has approximately 4.2 s between native-attempt duration and end-to-end
wall time, versus 1.1 s for A06. Split import/configuration/version/output
bookkeeping further before ascribing it to startup or promising its removal.

### 3. Heavy annotation: engine work cannot deliver the next large gain

In the strict-pass VEP WGS primary measurements, miss annotation consumes
**92.8% (A06) and 97.0% (C05)** of cached runtime. Making all non-annotation work
free would improve current cached runtime by only **1.077x and 1.031x**.
Halving it gives approximately 3.7% and 1.5% reductions in total runtime.

VCFcache can reduce the **effective time per original variant** by reusing more
work; ordinary routing improvements do not make the external tool annotate each
remaining miss faster. Larger gains here need fewer costly misses, amortized
annotator startup, or a separately optimized annotator recipe. Do not equate
cache hits with randomly sampled annotation cost: misses can be systematically
more expensive than cached variants.

Most promising follow-up: a small **cohort miss-union/batching** pilot with a
context-invariant, variant-local recipe. Annotate recurring misses once, preserve
each sample's own GT/FORMAT and metadata, and compare against an equally batched
uncached baseline. Count union/merge/output costs and startup amortization on
both sides. Cache-side batching is annotator-independent; no persistent service,
new cache format or embedded annotator is required. Context-sensitive VEP and
the failing vep-rs clinical recipe cannot simply be assumed safe for this.

An optional append-only, versioned cohort sidecar may capture recurrent misses
without changing the immutable base cache. Its gain depends on recurrence and
miss cost; the present campaign did not measure it. Do not claim that a hit rate
increase from 93% to 96% can only yield a few percent: it reduces the miss
fraction from 7% to 4%, potentially large when miss work dominates. Equally,
do not assume the remaining misses will recur or have equal costs.

Field-level/composable caches might eventually reuse expensive database payloads
across recipes without rebuilding every field, but would need dependency-aware
recipe signatures and exact field ownership. This is an interface/provenance
and correctness redesign, not an evidenced shortcut to lower latency here.
Defer it rather than expanding the next sprint.

## What the controls rule out or qualify

- **More RAM alone:** 8 GiB resident budget on core WGS still fell back after
  an 18.53-s native attempt; total 62.21 s versus 43.80 s adaptive. A late
  fallback is measurable waste. Better payload-aware preflight is low-risk;
  arbitrary larger buffers are not a scaling design. The buffer budget is not
  a process/cgroup memory limit.
- **Linear core scaling:** core WGS adaptive falls from 202.9 s at one CPU to
  62.4 s at eight and 43.8 s at 16: 4.63x, not 16x. FastVEP native WES is
  54.4 / 42.6 / 40.2 s at 1 / 8 / 16 CPUs. These are cached-pipeline controls;
  the direct arms were not rerun at each CPU count, so they are not matched
  per-core speedup comparisons.
- **Worker tuning:** fastVEP adaptive WES is 48.35 / 45.19 / 48.50 s at 2 / 4 /
  8 workers. This suggests a narrow heuristic test, not a universal setting.
  Native-request WGS worker cells actually fell back with adaptive's automatic
  worker selection; they do not test the requested native worker counts.
- **Skipping splitting:** fastVEP native WES changes from 40.23 to 38.04 s in
  single runs. Worth testing a verified already-split path, not dropping
  normalization safeguards based on a filename. VCF.gz is not consistently
  better than BCF. No case for a new on-disk format is established.
- **Instrumentation/noise:** three native fastVEP WES on/off pairs differ by
  -2.62%, +2.24%, +0.57%. This bounds neither all recipes nor adaptive overhead.
  A06 adaptive WES primary is 54.21 s versus 47.30 s uninstrumented; time-order,
  access state and instrumentation were not isolated for that difference.
  Use the uninstrumented references for the headline table; repeat only narrow
  promising controls before committing to an optimization.
- **CPU attribution:** perf was denied/unavailable. Coarse timing identifies
  useful stages but does not prove whether compression, allocation, memory
  bandwidth or actual storage reads dominate inside a stage. Measure before
  reaching for SIMD or custom low-level codecs.

## Do we need a deeper refactor?

**Not to release a useful, convincing product.** A whole rewrite, Parquet/new
key-value store or another annotator-specific backend is not justified by this
evidence. Large VEP gains already exist, fastVEP WGS is strong, and default
resident WGS is currently a fallback rather than an architectural advantage.

If the explicit target becomes **2–3x even for minimal vep-rs WGS**, a bounded
native routing/assembly design is a plausible next-stage prototype: process
coordinate-aligned chunks with bounded queues, avoid retaining the entire
annotated WGS, keep ordered output, and avoid repeating costly decode/overlay
passes. Still annotate the combined miss stream once; invoking an annotator
per chunk would multiply startup and can change context-dependent results.
BCF/CSI remains the external format. Preserve conservative fallbacks and strict
typed-record tests. This is materially more code and risk than heuristic tuning;
the expected benefit is not yet measured.

## Suggested decision gates and effort budgets

These are planning budgets, not completion promises.

| Priority | Work | Budget | Proceed only if |
| --- | --- | --- | --- |
| 1 | Audit/fix or replace the pinned rich vep-rs recipe; measure real information content | 1–2 engineering days for a decision | Complete-record parity; genuine useful extra annotations; no deliberately slow baseline |
| 2 | Payload-aware fallback prediction, focused window/worker/startup controls | 2–3 days | At least 5% reproducible end-to-end gain or removal of demonstrated fallback waste; no regressions |
| 3 | Core WGS assembly/index/copy prototype | 2–4 days | Move materially toward 2x on both samples at equal output/CPU/memory contracts |
| 4 | Small cohort miss-union experiment | 2–4 days, independent follow-up | Beat a fairly batched uncached baseline with strict per-sample equivalence |
| Deferred | General bounded-memory native redesign | 1–2 weeks to assess a prototype, not production-ready | Earlier gates identify enough removable cost; memory and correctness remain bounded |

Prefer one **3–5-day go/no-go sprint** covering priorities 1–2 plus one small
assembly experiment, not all rows simultaneously. Stop if gains are only a few
percent or require fragile special cases. The alternative of freezing now is
reasonable if the paper emphasizes annotation-rich pipelines and honestly shows
the fast-core break-even boundary; do not block release indefinitely on core WES.

## Before the full publication rerun

1. Freeze the engine/default/availability policy, source commit, helper build,
   recipes and all annotation resources. Release native as opt-in unless broader
   evidence justifies selection; document fallbacks and memory budgets.
   Ship a supported prebuilt helper/image if adaptive becomes the recommended
   or automatic path; otherwise ordinary installations can still execute legacy
   and never receive the advertised gain. Do not silently promise a fast path
   that requires an undisclosed optional build.
2. Restore an output-equivalent rich vep-rs condition or exclude that headline.
   Keep original VEP WES differences explicit; successful recovery WGS does not
   establish universal context independence.
3. Use ordinary real annotations for the main richness ladder; synthetic waits
   explain the model only. Record added INFO/CSQ fields, payload size and database
   operations rather than counting arbitrary plugin names.
4. Compare each recipe to its own efficient direct invocation with the same
   sorted/indexable output, cores and resource versions. Current direct controls
   use `vcfcache --uncached`, so include a minimal bare-annotator-plus-output
   control before presenting tiny-latency claims; shared wrapper overhead can
   affect ratios. Avoid untuned or OOM-constrained straw-man baselines.
5. Retain representative independent samples and a spread of hit rates, including
   the ~79% WGS case; do not advertise A06's 92.7% as universal. Honor the existing
   manuscript rule excluding gnomAD-contributing samples from gnomAD-cache
   performance figures, and audit provenance for final selected cohorts.
6. One run per independent sample/condition; small paired controls only for
   timing stability. Show per-sample distributions and cohort-stratified
   uncertainty, not a confidence interval from these two engineering genomes.
7. Include immutable cache build costs and cohort amortization; distinguish a
   bundled annotation cache from a locally built cache based on a bundled
   blueprint. The present timings exclude construction.
8. Capture each condition in its own Slurm step and choose wall limits from
   measured runtimes; keep long direct references separate or checkpoint-safe.
   Preserve previous publication figures and generate a new final output folder.

## Reproduction and evidence

Run from the repository root:

```sh
.venv/bin/python -m benchmarks.engine_profile.consolidate \
  benchmarks/engine_profile/consolidated-20261027
```

`cells.tsv` retains 94 original/recovery conditions separately; `attempts.tsv`
retains 95 actual state artifacts, including the superseded timeout's raw
`running` record (not a currently running job). `canonical_conditions.tsv`
selects the eight recovery configurations explicitly for the 86-condition view.
`matched.tsv`, `phases.tsv`, `instrumentation.tsv` and `targets.tsv` reproduce
the calculations. `SOURCE_MANIFEST.json` hashes the retrieved evidence.
Large variant outputs, input assets, perf binaries and per-second resource
streams remain remotely archived and are not duplicated here; final resources,
audit outcomes, commands, statistics and traces are included.

Earlier evidence: [branch evaluation](../../branch_eval/RESULTS.md),
[integration review and minimal clinical reproducer](../../branch_eval/integration-review-20261023/REVIEW.md).
The integration review's job-1328 status is an old snapshot, not a live status.

Consolidation validation: all seven TSVs and the source manifest reproduce
byte-for-byte; 86 canonical conditions, all 24 primary instrumentation-equivalence
checks, all eight recovery direct comparisons and actual WGS fallback labels
were checked programmatically. Repository tests: **563 passed, 21 skipped**;
Python lint passed. No new performance runs were submitted for this assessment.
