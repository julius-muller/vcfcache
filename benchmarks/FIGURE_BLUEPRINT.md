# User-facing figure blueprint

## Communication goal

This is a product benchmark, not a population-genetics study. The first figure
must let a laboratory or company answer this within seconds:

> My pipeline currently takes **X**, my expected cache hit rate is **f**, and I
> process **N** samples. How much waiting and compute will VCFcache remove?

The headline is **same annotations, less waiting**. Population distributions,
node metrics, confidence intervals, and cohort provenance remain available as
supporting evidence but do not lead the story.

## Hero figure: find your situation

### 1. One rule

Show the rule in plain language and small mathematical text:

`new time per sample ~= lookup time + (1 - hit rate) * current time`

Beside it, draw one full `current time` bar and four `time left with VCFcache`
bars. The filled part is time still spent; a pale segment with a large arrow is
time returned to the user.

| Typical situation | Hit rate | Approximate annotation time left* | Approximate speedup* |
|---|---:|---:|---:|
| Distantly related WGS | 50% | 50% | 2x |
| Related/cohort WGS | 80% | 20% | 5x |
| WES | 90% | 10% | 10x |
| Very high reuse | 95% | 5% | 20x |

`* plus the small measured lookup and preprocessing overhead`

The labels are examples, not biological guarantees. The measured external WGS
and assay distributions determine the final ranges displayed around them.

### 2. Fast pipeline or slow pipeline

Use a three-column lookup strip for one selected hit rate. For an 80% hit rate:

| Your pipeline today | Approximate time with VCFcache* | Time returned per sample* |
|---:|---:|---:|
| 10 minutes | 2 minutes | 8 minutes |
| 1 hour | 12 minutes | 48 minutes |
| 10 hours | 2 hours | 8 hours |

This makes the distinction intuitive: relative speedup is driven mostly by hit
rate, while absolute compute savings become larger for expensive pipelines.

### 3. Ten samples or one thousand

Show three sample-count cards: `10`, `100`, and `1,000`. Each card reports:

- original total runtime;
- runtime with VCFcache;
- total time or compute returned;
- effective per-new-sample runtime.

For a bundled Zenodo cache the user has no cache-build cost, so the benefit
starts with the first sample. For a cohort-built cache, show one additional
break-even marker calculated from the measured build time. Do not imply that
per-sample speedup grows with cohort size: cohort size multiplies total savings
and amortizes only the one-time build cost.

### 4. Trust badge

Place a visually separate check-mark card beside the performance graphic:

> **Same annotations**
> Every cached result was compared with its uncached result. Any annotation
> difference fails the task.

The small print documents the known VEP 115.2 `HGNC_ID` nondeterminism from
Ensembl issue 1959. This is an upstream VEP behavior observed between VEP runs,
not an annotation change introduced by VCFcache.

## Supporting evidence

These plots substantiate the simple hero graphic but should appear after it:

1. **Real assays:** individual relative runtimes for 20 panel, 20 WES, and 49
   completed WGS pairs using the bundled Zenodo AF >= 1% cache.
2. **Independent WGS:** all 52 KPGP, SGDP, and PGP evaluation genomes comparing
   bundled AF >= 10%, bundled AF >= 1%, and disjoint three-genome cohort caches.
   "Evaluation" means that none of these genomes was among the three genomes
   used to construct its cohort cache; it does not assert individual-level
   exclusion from the independently published gnomAD source population.
3. **Pipeline cost:** one representative sample run with vanilla VEP, VEP
   `--everything`, and two calibrated `SyntheticDelay` settings.
4. **Observed versus predicted:** measured relative runtime against hit rate,
   with the simple runtime rule overlaid. This demonstrates that users can apply
   the calculator to their own runtime.

## Rendering rules

- Produce every publication and repository plot in R with ggplot2. Python may
  collect and validate source tables, but it is not the plotting backend.
- Lead with minutes/hours saved, not p-values or hardware counters.
- Use one consistent color: dark for time still spent, light green for time
  returned, and a check-mark color for correctness.
- Always show both relative runtime and a concrete time example.
- Label modeled values as estimates and measured values as observations.
- Never pool assembly-specific cohorts into one unexplained average.
- Keep raw points, exact commands, resource metrics, and provenance in the
  downloadable benchmark package.

## Execution design

Use one validated paired or multi-condition execution per sample, with balanced
condition order. The biological/input sample is the independent unit. The
existing HG02079 repeat differed by 0.84% uncached and 0.08% cached, so complete
cohort warm-ups and three technical reruns were removed before collection.

Do not resize workers during the primary/external comparison. A controlled
pipeline-cost campaign may use different hardware only when all its scenarios
run on that same hardware and are presented as a separate experiment.
