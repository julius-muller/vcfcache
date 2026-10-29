# Final performance decision analysis

29 October 2026 · frozen campaign evidence · no additional annotation submitted

## 1. Decision in one paragraph

**Proceed with the current adaptive implementation as the release candidate for
annotation-rich workflows and WGS; do not delay the paper for a whole-engine
rewrite.** The complete campaign supports substantial WGS acceleration with VEP,
SnpEff and fastVEP, including a small but consistent WGS benefit with core vep-rs
in the matched extension. It does **not** support universal acceleration of fast
annotators at WES. Keep native available as an opt-in, measured specialist path,
not as an abandoned experiment or an unconditional default. Further work should
be a separately approved, bounded routing/assembly experiment, or cohort-level
miss batching with an equally batched direct control—not another broad tuning
matrix. Clinical vep-rs remains timing-only and must not supply a validated
headline. The completed manuscript campaign does not need another full rerun
merely to reach this decision.

Read [DECISIONS.md](DECISIONS.md) for concrete next actions and gates. The
[figure gallery](README.md#figures) and [captions](FIGURE_LEGENDS.md) accompany
this analysis. Every number below is recoverable from the linked TSVs; no
technical repeat averaging or fastest-run selection was used.

## 2. What actually finished

All annotation work and final collectors finished on **29 October at 12:57 UTC
/ 13:57 CET**, before the Monday 2 November target. A fresh queue inspection
found no queued/running jobs and all six benchmark workers idle. That idle state
is intentional: the approved matrix is finished, not waiting for a dispatcher.
Post-processing is local active work, not another Slurm benchmark.

| Inventory | Number | Interpretation |
|---|---:|---|
| Canonical terminal conditions | 600 | Distinct planned sample/recipe/cache conditions |
| Archived execution attempts | 604 | Four original header failures retained, then superseded |
| Plotting views | 612 | Includes 12 aliases of existing measurements, not extra observations |
| Distinct genome identifiers overall | 54 | Not 600 biological samples |
| Accepted under declared comparators | 572 | Several comparison scopes; not 572 fresh strict direct comparisons |
| Clinical vep-rs timing-only | 24 | 12 WES and 12 WGS conditions with CSQ differences |
| SnpEff strict content failures | 4 | Known warning discrepancy; excluded from accepted speedups |

The INFO-header repair completed all **36** affected GRCh37 clinical SnpEff
slots: 34 passed and two retained the known content discrepancy. Across core and
clinical SnpEff together, four strict discrepancies remain; they are not ongoing
jobs or unexplained process crashes. All 600 conditions have successful cached
timing measurements, but successful execution alone does not establish output
equivalence.

| Experiment family | Canonical conditions | Distinct genomes | Accepted | Purpose |
|---|---:|---:|---:|---|
| Full WGS cache strategies | 300, plus 12 plot aliases | 52 | 300 | VEP and basic fastVEP; three cache strategies |
| Clinical assay series | 180 | 24 | 176 | fastVEP and SnpEff; panel/WES/WGS |
| Matched VEP assay series | 36 | 12 | 36 | Common GRCh38 assay comparisons |
| Synthetic pipeline spectrum | 12 | 1 | 12, semantic/composed scope | Six plugin waits, two cache hit rates |
| Pinned vep-rs extension | 72 | 12 | 48 | Core and clinical; three assays |

Families overlap in donors and must not be added as independent samples. The
52-genome strategy cohort comprises **20 KPGP, 20 SGDP and 12 PGP** genomes.
The common GRCh38 assay comparison comprises **six KPGP and six SGDP** genomes.
The separate PGP assay series has 12 GRCh37 genomes. There is no 1000 Genomes
test cohort in this final primary matrix.

### Execution and comparison controls

- Primary requested **and effective** engine: **adaptive**, eight allocated and
  affinity-bound CPUs, node-exclusive execution, 28-GB job allocation. CPU
  allocation does not imply all eight cores stayed busy.
- 562 conditions used `fe197f34d`, 36 the header repair `787a3081a`, and two
  carried-forward checks `b4b43512c` with identical frozen product-file hashes.
  Helper versions are retained per row, not silently conflated.
- **252** conditions use newly measured paired direct controls: 248 within their
  final attempt and four successful original direct controls reused after the
  header repair. The other **348** retain compatible historical direct timings.
  Reuse avoids unnecessary annotation but does not create a randomized,
  contemporaneous direct-versus-cached trial.
- Each fresh direct/cached timing has a separate Slurm step/cgroup. Aggregate
  CPU counters and condition-specific cgroup peak memory are available for those
  controls. Historical denominators are not given invented fresh resource data.
- Recipes, parameters, input/cache hashes, annotation assets, commands and
  comparator decisions are frozen. The inventory records VEP 115.2, fastVEP
  0.3.0, SnpEff 5.4 and the pinned vep-rs recipe requirement 116.2; the binary
  fingerprints, not a claim about current upstream releases, define the runs.
- Panel and WES conditions are **interval-restricted WGS callsets**, not separate
  wet-lab panel/WES sequencing cohorts. Different assemblies, provider pipelines,
  annotation database contents and storage paths are not interchangeable.

Source: [canonical cells](tables/cells.tsv), [attempts](tables/attempts.tsv),
[audit](tables/FINAL_AUDIT.json), [cache/recipe inventory](tables/cache_inventory.tsv).

## 3. Strongest common-cohort result: workload-specific acceleration

The following medians use the same 12 GRCh38 donor identifiers for each assay,
the gnomAD AF ≥1% cache blueprint, and each recipe's **own** direct baseline.
Clinical fastVEP and SnpEff add population/clinical annotations; they are not
information-equivalent implementations of VEP `--everything` or core vep-rs.
The table therefore compares the benefit of caching, not raw annotator rankings.

| Recipe | Panel speedup | WES speedup | WGS speedup | WGS samples faster |
|---|---:|---:|---:|---:|
| VEP everything | 2.59× | 2.46× | **10.23×** | 12/12 |
| fastVEP clinical | 1.00× | **0.81×** | **2.88×** | 12/12 |
| SnpEff clinical | 1.27× | 2.65× | **10.81×** | 12/12 |
| vep-rs core | 1.33× | **0.46×** | **1.32×** | 12/12 |

| WGS recipe | Median speedup, 95% bootstrap interval | Median paired time saved | Median runtime reduction |
|---|---:|---:|---:|
| VEP everything | 10.23× [8.73, 10.38] | **2.30 hours/sample** | 90.2% |
| fastVEP clinical | 2.88× [2.78, 2.95] | **5.11 minutes/sample** | 65.3% |
| SnpEff clinical | 10.81× [9.81, 11.59] | **2.34 hours/sample** | 90.7% |
| vep-rs core | 1.32× [1.22, 1.37] | **16.4 seconds/sample** | 24.2% |

This is a credible user-facing story: reuse avoids a large amount of repeated
work in annotation-rich pipelines, and still helps WGS with a fast annotator.
It is not honest to call the core-vep-rs WGS improvement “stunning,” or to hide
its WES loss behind the much larger VEP result.

### The losses are important product information

Clinical fastVEP WES is faster with caching for only **2/12** genomes, with a
median **7.48-second penalty**; the median speedup interval is [0.78, 0.94].
Core vep-rs WES loses for **all 12**, with a median **11.47-second penalty** and
speedup interval [0.44, 0.55]. Panel clinical fastVEP is approximately break-even
(6/12 faster), not a robust win. VEP Panel has one slower sample despite its
strong median. Do not drop these points or promise “always faster.”

The apparently better core-vep-rs Panel ratio does not contradict the WES loss:
the direct tool's startup and small-input execution behavior are not a simple
linear extrapolation from WES, while the cache engine switches sparse access
strategies. A ratio alone does not establish the startup mechanism.

The 95% intervals resample **donors within cohort**, preserving the six/six
allocation. They describe sample-to-sample uncertainty conditional on these
cohorts/settings, not repeated-run timing variability or representativeness of
all laboratories. One timing per condition means order, node and OS-cache effects
remain potential influences. Closely related donors, if present, also weaken a
strict independence interpretation; these are descriptive cohort estimates,
not a population-genetics inference.

Source: [matched summary](tables/primary_grch38_summary.tsv),
[cohort-stratified tables](tables/assay_summary_by_cohort.tsv),
[lowest-speedup samples](tables/lowest_speedup_samples.tsv), Figure 1.

## 4. All 52 WGS genomes: cache strategy matters

This family uses VEP everything and **basic fastVEP**, not the clinical fastVEP
recipe in Section 3. “AF ≥1%” and “AF ≥10%” describe variant selection in the
gnomAD blueprint, not a target sample's observed hit rate.

| Cohort | Cache | Median hit rate | VEP speedup | Basic fastVEP speedup |
|---|---|---:|---:|---:|
| KPGP, n=20 | gnomAD AF ≥10% | 81.3% | 4.84× | 3.00× |
| KPGP | gnomAD AF ≥1% | 92.0% | **8.82×** | **3.57×** |
| KPGP | Three held-out cache seeds | 82.9% | 4.99× | 2.89× |
| SGDP, n=20 | gnomAD AF ≥10% | 83.3% | 5.63× | 2.53× |
| SGDP | gnomAD AF ≥1% | 93.0% | **10.76×** | **3.50×** |
| SGDP | Three held-out cache seeds | 81.5% | 4.80× | 2.47× |
| PGP, n=12 | gnomAD AF ≥10% | 84.0% | 6.59× | 2.81× |
| PGP | gnomAD AF ≥1% | 92.7% | **13.15×** | **3.61×** |
| PGP | Three held-out cache seeds | 77.5% | 4.71× | 2.72× |

Paired AF ≥1% versus AF ≥10% cached-runtime improvements are approximately
**1.98× / 2.02× / 1.65× for VEP** and **1.24× / 1.44× / 1.23× for fastVEP**
in KPGP/SGDP/PGP. These are medians of within-donor ratios, not ratios of the
medians above. AF ≥1% improves every VEP sample here, but only 83–90% of basic
fastVEP samples; higher hit rate does not guarantee faster total runtime once
lookup/output costs and run-to-run effects are included.

The three-genome cache is useful but is **not consistently better** than the
AF ≥10% cache: it is broadly similar for KPGP and generally less favorable for
SGDP/PGP in this matrix. Do not select a cohort cache just because its name
sounds more personalized. Consider its actual hit rate, build cost and recipe.

### Realistic low-hit cases are retained

SGDP **LP6005592-DNA_C05** has approximately **60.9%** hits against AF ≥10% and
**57.9%** against the three-genome cache. Basic fastVEP speedups are **0.935×**
and **0.970×**, respectively: caching loses slightly. The same genome with the
AF ≥1% clinical configuration has **79.0%** hits and **2.21×** fastVEP speedup;
core vep-rs reaches **1.095×**. These are different recipe/cache conditions,
not evidence that ancestry alone explains a change in speed.

Thus, neither “every WGS is 90% hits” nor “caching always wins at 60% hits” is
supported. Blueprint breadth, variant representation/filtering, provider
pipeline, assembly, annotation cost and hit composition can all matter. This
campaign does not isolate ethnicity, relatedness or sequencing technology as
causal explanations. Jaccard overlap between two samples is also not the same
quantity as a sample's hit fraction against a union/population cache.

### What “held-out” does and does not mean

All three cohort-cache seed sets are ID-disjoint from their own 52-genome-family
test donors; [the explicit seed audit](tables/cohort_seed_audit.tsv) records the
nine seed identifiers and test counts. This is a held-out **cache-training**
test. It does not prove absence from gnomAD's donor set or that all donors are
unrelated. No gnomAD donor-disjointness claim is made. These are external public
cohorts, not a 1000 Genomes benchmark presented as independent of gnomAD.

Source: [WGS summaries including confidence intervals](tables/wgs_strategy_summary.tsv),
[paired strategy contrasts](tables/paired_strategy_contrasts.tsv), Figure 2.

## 5. Output fidelity: a necessary limit on the claims

“572 accepted” is useful bookkeeping, but the paper must describe the actual
comparison scopes:

| Scope | Conditions | Accepted interpretation |
|---|---:|---|
| Retained full cached reference plus annotator-specific direct comparison | 336 | 156 fastVEP strict direct comparisons; 180 VEP direct semantic comparisons with the predeclared HGNC_ID handling |
| Fresh paired complete-record/header comparison | 228 | 224 passed; four SnpEff warning differences failed |
| Clinical vep-rs paired output with CSQ differences | 24 | Timing only, not validated acceleration |
| Composed vanilla-VEP semantic reference for delayed pipelines | 12 | Semantic/composed assurance, not per-delay retained full-record direct validation |

Complete-record comparison retains INFO/FORMAT values and relevant headers;
only declared ordering canonicalizations are permitted. The generic fastVEP
comparator's name does not mean that only CSQ was checked for SnpEff. No new
HGNC_ID exception was introduced for fastVEP or vep-rs. Full BCF outputs remain
on the cluster; this compact local archive contains commands, comparison
results/examples, checkpoints, metrics, recipes and hashes, not the full BCFs.

### Clinical vep-rs

The 12 WES conditions have **5–34** differing records each; WGS has **13–126**.
Relevant headers and the non-CSQ comparison pass. Median timing ratios are
**1.93× WES** and **3.18× WGS**, but neither belongs in the validated performance
headline. Earlier multi-engine checks and a cache-independent reproduction
support an annotator/subset-context issue rather than a routing regression.
They do not establish that the emitted annotation is interchangeable.

The pinned annotator is intentionally not patched/upgraded during this sprint.
Keep these timings visible in a diagnostic-only figure and exclusions table.
An independent future recipe investigation can be proposed if a validated rich
vep-rs headline is indispensable; do not silently “fix” the comparator.

### SnpEff

The four failures are **hu6664F3**, core/clinical × WES/WGS, each one record at
**chr2:11273459 C>T**. Saved examples identify the discrepancy as
`WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS` in ANN transcript entries. The known
historical failure status is retained. All four remain strict failures and
are excluded from accepted speedup summaries. The compact examples are
truncated; a full-record reinspection needs the retained remote BCFs.

### VEP

The existing, explicitly documented HGNC_ID semantic exception is retained
where already approved; the comparator also canonicalizes ordering within the
CSQ DOMAINS field. The audit retains those rules and whether an approved cached
semantic exception was needed. It exports raw direct differences rather
than implying byte-for-byte identity: 96 canonical VEP conditions contain raw
differences across the full families, with 14,036 **condition-record comparison
differences**, not that many unique biological variants. Thirteen matched VEP
assay conditions account for 148 such differences. The older manuscript's
“265 differences” belongs to its original comparison set; it is not a count
that can simply be pasted onto this expanded final campaign. Preserve both
sets' provenance and explain the current comparison scope in revised methods.

The appropriate claim is **reuses annotations under a frozen, compatible recipe,
with equivalence evaluated under declared comparisons**, not an unconditional
promise that all annotators emit identical output on arbitrary subsets.

Source: [per-condition comparator audit](tables/comparator_audit.tsv),
[excluded rows](tables/validation_exclusions.tsv), Diagnostic D1.

## 6. Measured bottlenecks: where further work can pay

These are medians of actual phase wall times, not inferred intercepts from two
input sizes. Per-phase medians do not necessarily add to the median total. The
per-condition decomposition does add exactly; parallel worker durations are not
summed into fictitious elapsed time. “Annotate misses” includes invoking the
external recipe and its I/O/startup, not only pure consequence calculation.

| GRCh38 WES recipe | CLI/other | Routing | Assembly | Annotate misses | Total cached |
|---|---:|---:|---:|---:|---:|
| VEP everything | 1.20 s | 19.23 s | 3.42 s | 135.61 s | 164.31 s |
| fastVEP clinical | 1.10 s | 16.47 s | 2.50 s | 37.11 s | 59.70 s |
| SnpEff clinical | 1.19 s | 13.94 s | 1.74 s | 180.44 s | 199.37 s |
| vep-rs core | 1.12 s | 11.47 s | 1.34 s | 7.25 s | 22.49 s |

**The next WES lever is routing, not another generic Python-startup campaign.**
The large genome-wide fixed term that motivated adaptive restriction has already
been addressed. Remaining sparse-query/candidate reads, indexed access, decoded
payload handling, subprocess boundaries and output construction are candidates
for a new phase-controlled experiment—not proven removable costs.

### What would it take to beat uncached core vep-rs at WES?

Holding measured miss annotation unchanged, median cache-side work must fall
by **77.4%** merely to reach parity. Halving all cache-side work yields only
**0.71×**, and cutting it by 75% yields approximately **0.97×**. Removing the
entire measured route phase gives **0.98×**; removing both route and assembly
gives **1.10×**. Even eliminating every non-annotation cost gives only **1.55×**
median speedup with the present miss stage.

These are counterfactual bounds, not physical impossibility proofs: a redesigned
workflow may also change invocation cost. Nevertheless, a quick 5–10% tuning
patch cannot meet this goal, and “stunning core WES acceleration” is not a
credible near-term commitment on the current evidence.

For clinical fastVEP WES, the hurdle is smaller: **42.4%** median cache-side
reduction to parity. Halving that work predicts **1.05×**, while deleting all
non-annotation work gives **1.47×**. Core vep-rs and fastVEP therefore must not
be treated as the same optimization target.

### WGS has materially better fast-annotator headroom

| Recipe | Current WGS speedup | If half of cache-side time vanished | If all vanished |
|---|---:|---:|---:|
| VEP everything | 10.23× | 11.14× | 12.21× |
| fastVEP clinical | 2.88× | **3.89×** | 5.97× |
| SnpEff clinical | 10.81× | 11.29× | 11.81× |
| vep-rs core | 1.32× | **2.13×** | 4.77× |

The last two columns remove costs mathematically while holding miss work fixed;
they are not measurements or promises. They show where engineering could have
meaningful leverage: fastVEP WGS and core-vep-rs WGS, rather than already
miss-dominated heavy pipelines. Current WGS route/assembly medians are
**34.8/32.4 s** for clinical fastVEP and **20.1/18.2 s** for core vep-rs.

### Heavy pipelines: reduce expensive work, not only wrapper costs

Miss annotation occupies a median **84.7%** of cached VEP WGS and **91.7%** of
clinical SnpEff WGS. Deleting all non-annotation work improves their cached
runtime by only about **1.18× / 1.09×**, respectively. In the largest synthetic
wait cases the annotation share rises to about 98–99%.

Routing optimization does not reduce the annotator's intrinsic per-miss work.
It reduces effective time per *original input variant*. Larger gains need fewer
expensive misses, safely deduplicated repeated work, or amortized invocation
cost—not an unsupported claim that VCFcache accelerates every database lookup
inside the external annotator.

The simple model remains useful, but distinguish
`T_cached − (1−f) × T_direct` from measured overhead: this residual also contains
startup and systematic differences between cached and missing variant costs.
Effective microseconds per record are recipe/scale-dependent totals, not universal
`t_ann` constants. A two-point fit does not prove the mechanism of a fixed cost.

Source: [measured phases](tables/phase_seconds_median.tsv),
[individual-stage bounds](tables/individual_stage_what_if.tsv),
[headroom curves](tables/cache_side_what_if.tsv),
[target-speedup requirements](tables/target_speedup_headroom.tsv), Figure 3.

## 7. Pipeline complexity and absolute time saved

The synthetic spectrum uses one real WGS genome, **KPGP-00319**, at **80.23%**
and **90.26%** observed hits. The self-written plugin waits 0.5, 1, 2, 4, 7 or
10 milliseconds **per consequence**, not per variant. The inputs are two cache
strategies on the same genome, not two independent cohorts.

At the 10-ms setting the historical direct runtime is **22.71 hours**. Current
cached runtime is **4.49 hours** at 80.23% hits (**18.22 hours saved**, 5.06×),
and **2.24 hours** at 90.26% (**20.46 hours saved**, 10.13×). This is a clear
illustration of large absolute savings for costly annotation. At the 0.5-ms
setting the direct runtime is 4.21 hours and the AF ≥1% cached runtime 0.47 hours.

The points are not perfectly monotonic in speedup: direct controls are historical,
one timing per setting, and consequence counts/costs and machine state are not
captured by a one-parameter per-variant model. Do not fit a smooth claimed law
through these points or label simulated waiting as ten real database plugins.
The semantic/composed validation limit in Section 5 must accompany this figure.

These scenarios demonstrate delay-related time saving, **not** whether ten real
fastVEP databases become free by adding ten cores. That requires a separate,
matched database-count × CPU-allocation experiment with direct and cached arms.
The existing thread-limit pilots and extension recipes do not settle it.

Source: `synthetic` tibble in [manuscript_data.rds](tables/manuscript_data.rds), Figure S1.

## 8. Memory, CPU cost, setup and cohort amortization

### Faster does not always mean lower memory

The fresh matched GRCh38 WGS controls give median cached/direct aggregate CPU
ratios of **0.36** for clinical fastVEP, **0.12** for clinical SnpEff and **0.74**
for core vep-rs. Those are useful compute savings. But their median peak-memory
ratios are approximately **3.24×, 1.29× and 1.31×**, respectively. Absolute cached
peaks are approximately **11.6, 17.5 and 8.4 GiB**. Cgroup peak includes aggregate
charged memory, including file-backed/page-cache effects; it is not a single
process's RSS or proof that all of it is irreducible application heap.

Clinical fastVEP WES costs about **2.00×** direct CPU, and core vep-rs WES about
**2.25×**. These losses should inform an eventual bypass recommendation, not only
the runtime figure. The allocation was 28 GB; parallel-density recommendations
need measured peaks and contention tests, not an assumption that eight workers
scale linearly. `samples_per_day` is a serial-slot latency extrapolation, not a
measured multi-sample throughput campaign.

Source: [fresh paired resources](tables/paired_resources.tsv),
[resource summaries](tables/resource_summary.tsv), Figure S3.

### Cache origin is not interchangeable

VEP population caches are ready-made Zenodo assets. The basic fastVEP caches were
**locally annotated from the bundled gnomAD or cohort blueprint**, not downloaded
as already annotated fastVEP caches. Clinical recipes likewise require their
own recipe-compatible payloads; a matching genome/AF threshold is insufficient.
The [inventory](tables/cache_inventory.tsv) records actual recipe and payload
hashes and any available build provenance. Zenodo DOI/source aliases remain in
[the provenance/build table](tables/recorded_build_amortization.tsv).

The audit found that historical strategy metadata and the actual fastVEP payload
provenance describe slightly different build durations. We use the provenance
whose **cache SHA256 matches the executed payload**, not the earlier strategy
number or a near-zero reuse duration. Both numbers remain inspectable. This is a
provenance correction in this analysis, not a mutation of the historical files.

### What can be said about amortization?

Payload-bound historical build costs are about **3.4–10.4 min** for the three
basic fastVEP cohort caches and **1.85–3.95 h** for the VEP cohort caches. Basic
fastVEP AF ≥1% builds are about **8.7 min (GRCh37)** and **25.1 min (GRCh38)**.
Using each cohort's observed median paired savings, the recorded annotation build
is repaid after approximately **2–3 new WGS samples for fastVEP cohort caches**,
**4–6 for fastVEP AF ≥1%**, and **two for VEP cohort caches**.

These are **conditional historical-cost calculations**, not measured end-to-end
cache preparation on today's hardware. FastVEP build costs start from an available
blueprint; download and blueprint preparation are excluded. Ready-made VEP caches
need no local annotation build, but download/setup time is unmeasured, not zero.
The small repayment counts refer to *additional test samples after preparation*,
not a recommendation to build a cache for a two-sample lab.

Figure S5 explicitly assumes normalized build costs to illustrate amortization
without pretending they were measured. Figure S6 combines recorded payload-bound
historical fastVEP build costs with the current **mean** direct/cached sample
costs; it is a conditional additive cohort scenario, not a new cohort benchmark.
Per-condition medians, ratios of medians and aggregate throughput ratios are
different statistics and are not substituted for one another.

## 9. Main versus adaptive versus native: benefit and maintenance

This complete campaign directly evaluates **adaptive**, not all engines in a
new paired full-cohort trial. The earlier branch and engine profiles remain
relevant engineering evidence but must not be pooled into the new confidence
intervals: they used fewer genomes, different CPU settings and sometimes native
requests that actually fell back to adaptive.

| Path | Current assessment | Recommendation |
|---|---|---|
| Legacy/main behavior | Established compatibility; substantial old small-input fixed cost | Keep as compatibility fallback/reference, not the performance headline |
| Original fused | Dense-input benefit but severe single-contig/small-input pathology | Do not select as a universal default |
| Adaptive | Broad applicability, explicit sparse/dense restriction, complete 600-condition evidence | Production candidate and manuscript primary |
| Resident native | Real measured WES benefit in earlier controls; more state/memory and an additional implementation | Retain opt-in, with parity tests, memory guard and actual-backend reporting |

Native's WES improvement is worth preserving. It is not, however, evidence that
every WES workload beats direct annotation, or that a WGS result labeled
“requested native” ran natively. Earlier pilots include native fallbacks; use the
[corrected backend audit](../../branch_eval/integration-review-20261023/cells_backend_audit.tsv).
The existing implementation shares orchestration but still adds a resident Rust
path and Python adapter, with their own header, allele-pairing and memory
invariants. Its incremental gain needs to pay for those tests and maintenance.

The 36-cell incremental sprint is already complete. Indexed spooling, smaller
native query windows and trusted-hit assembly did not show consistent total-time
gains. Bounded native routing improved sparse fastVEP WES by roughly 5–8% but
regressed dense core-vep-rs WGS by 23–31%. Early memory preflight avoided about
22–30% wasted fallback time in the deliberately over-budget native controls;
it does **not** accelerate an already adaptive condition. Those observations do
not justify enabling all experimental switches to advertise a compounded gain.

There is a release-interface gap: **the CLI still defaults to `legacy`**;
`--engine auto` selects available adaptive support or falls back to legacy, not
native. The measured results must be accompanied by `--engine adaptive` or a
reviewed default/packaging change before ordinary users can be promised the
same path. This analysis deliberately changes neither defaults nor packaging.

Source: [screening conclusions](../SCREENING_RESULTS.md),
[previous engineering assessment](../../engine_profile/consolidated-20261027/ASSESSMENT.md),
Figure S4 (historical comparison, not fresh causal A/B evidence).

## 10. Recommended next decisions, ranked by value

1. **Use these results and freeze the manuscript dataset.** Lead with WGS reuse
   across recipe complexity, retain the unfavorable fast WES cases, and describe
   semantic scopes precisely. No further full campaign is needed unless product
   behavior or recipe/output semantics materially change.
2. **Review adaptive availability/defaults and native support as a release task.**
   Package the compatible helper; test auto/fallback/no-helper paths, receiving
   headers, multiallelics, duplicate loci and actual-backend reporting. Keep
   native opt-in. A small representative release check is more proportionate
   than another 600-condition run for an interface-only change.
3. **If more speed is desired, approve one bounded WGS routing/assembly pilot.**
   Use the ordinary and low-hit SGDP genomes, core vep-rs and clinical fastVEP,
   frozen recipes/output formats, matched CPUs and alternating control order.
   Target >=10% reproducible total-time improvement without >5% reproducible
   regression, strict record/header parity and an explicit peak-memory cap.
   A 2× core-vep-rs WGS result requires a large reduction; it is a target, not
   an already demonstrated result. Do not rerun the full cohort before this gate.
4. **For larger heavy-pipeline gains, test cohort miss batching separately.**
   Annotate the union of recurring misses once, redistribute only safe
   variant-local annotation while preserving per-sample INFO/FORMAT/GT, and
   include union/merge/output costs. Compare against an equally batched direct
   control. Start with a recipe whose subset invariance passes; known clinical
   vep-rs/SnpEff/VEP context differences cannot simply be assumed safe. No embedded
   annotator, persistent service or new storage format is required for a pilot.

A self-learning sidecar can be a later append-only, versioned layer, leaving
the DOI-pinned base immutable. Its value depends on **cost-weighted miss recurrence**,
not hit rate alone. A change from 93% to 96% hits cuts the miss count from 7% to
4%; when miss work dominates, that can matter greatly. Conversely, it helps
little if the expensive misses are private or the annotator startup dominates.
Neither recurrence nor the resulting gain was measured here.

Do **not** start a new key-value store, Parquet/Manticore conversion, full Rust
rewrite or annotator integration merely because another implementation is fast.
The current evidence identifies stages, not a demonstrated storage-engine
solution. Those changes add format/compatibility burden and are outside the
agreed sprint scope. Reconsider a deeper refactor only after a bounded prototype
shows a meaningful end-to-end gain and an affordable correctness surface.

**Bottom line:** there is already a publishable, useful product story, and there
is still targeted fast-annotator WGS headroom. There is no evidence-based reason
to postpone everything until consequence-only WES becomes a large win.

## 11. Reproducibility and caveats to carry into the paper

- Frozen source: [SNAPSHOT.json](source/SNAPSHOT.json); its archive SHA256 is
  `9279773a1beac394fa404c6855002aaa4cf956656bd87fd2fff843edd3b9df74`.
- Archive and per-condition hashes allow re-auditing commands, recipes, parameter
  files, comparison scope, source build and metrics without retrieving large BCFs.
  Full-record reinspection still requires the retained remote outputs.
- [The R entry point](../analyse_final.R) creates all figures using ggplot2;
  [the helper](../final_analysis_helpers.R) supplies cohort-stratified bootstraps
  and file-preserving exports. The named tibble collection is
  [manuscript_data.rds](tables/manuscript_data.rds).
- Bootstrap intervals are conditional and descriptive; no multiple-comparison
  significance claims or timing-variance estimates are made. Related donors and
  systematic historical-versus-current execution differences are not removed by
  resampling. Synthetic scenarios have one donor and no inferential interval.
- Fixed eight-CPU budgets do not establish linear threading or cluster throughput
  scaling. Earlier configured-thread pilots are not strict CPU-quota experiments.
- Historical-to-current cache gains are useful operational context, not a clean
  causal estimate of a single code change. Different storage placements also
  prevent attributing cross-family differences solely to engine/annotator design.
- No unsupported “zero result changes,” “all WGS hit rates above 90%,” “all donors
  outside gnomAD,” “all setup is negligible,” or “all fast annotators benefit at
  every assay size” statement should enter the revised manuscript.
