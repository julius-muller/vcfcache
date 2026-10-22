# vep-rs in the annotation-heavy regime — genome-wide gnomAD

Slurm job 1251, 2026-10-22, sl-w1, exclusive, 16 vCPU. 66 cells, all `rc=0`,
every level's defining field asserted present in its output.

Job 1242 measured vep-rs on transcripts alone and found it far below break-even.
That was its cheapest possible state. This is the other end.

## The fits

Two input sizes, three reps each, median per size, `t = s + c·N`.
Thresholds: **29.7 µs/variant** legacy, **15.0 µs/variant** fused-routing
(`../fused-routing-threshold-check/`).

| rung | per variant | rec/s | legacy | fused |
|---|---|---|---|---|
| core only | 7.42 µs | 122,629 | lose | lose |
| + genome-wide gnomAD AF | **31.95 µs** | 30,674 | **win** | **win** |
| + ClinVar (= clinical level) | 40.76 µs | 24,378 | win | win |
| + HGVS/symbol/canonical | 52.93 µs | 18,963 | win | win |
| + transcript scores | 76.36 µs | 13,243 | win | win |
| `--everything` + clinical | 59.65 µs | 16,880 | win | win |
| `--everything` + heavy | 95.59 µs | 10,559 | win | win |

**One genome-wide supplementary source is enough.** It moves vep-rs from 7.42 to
31.95 µs/variant — past even the legacy threshold, without touching the fused
engine. My earlier "two to three sources" estimate was wrong, and wrong for an
instructive reason: it came from the exome-restricted file, which has almost no
allele frequency to return outside exons.

Only bare core stays below. That matters: the rule keeps a real negative region
and so remains falsifiable.

## Genome-wide costs less than exome-restricted and delivers six times more

| source | records | per variant | rows annotated (of 3,787,554) |
|---|---|---|---|
| genome-wide joint | 99,942,855 | 31.95 µs | **3,540,863** (93.5%) |
| exome-restricted | 183,717,261 | 48.92 µs | 573,173 (15.1%) |

The exome file is both **53% more expensive per variant and returns 6.2× fewer
annotations**. It holds nearly twice as many records packed into exonic clusters,
so each tabix query decompresses a larger block, and off-target queries pay that
cost to return nothing. Using it would have produced a slower, emptier, entirely
defensible-looking clinical arm.

## Two caveats that bound how precisely these numbers can be quoted

**The source-count sweep is sub-additive, so it is a lower bound.** Marginal cost
per instance falls with instance count:

| instances | total | per source |
|---|---|---|
| 1 | 31.95 µs | 24.53 µs |
| 2 | 53.63 µs | 23.11 µs |
| 4 | 68.56 µs | 15.29 µs |
| 8 | 98.66 µs | 11.41 µs |

All instances point at the *same* file, so later ones hit a warm page cache. A
real recipe over gnomAD **plus** SpliceAI **plus** CADD **plus** REVEL cannot
share cache that way and would cost more than this sweep implies. The sweep
models the floor of multi-source cost, not its expectation.

**Four rungs fit a negative start-up** (−0.060 to −0.278 s), which is physically
impossible and means `s + c·N` is the wrong shape for them. The small slice is
1.5 Mb of dense lookups; the big slice spans all of chr22, so per-variant cost
is not constant across the two. For those rungs the two-point fit *overestimates*
c. Every one of them sits far enough above both thresholds that the verdict is
unaffected, but the individual µs figures should be read as upper bounds rather
than quoted as measurements. A third input size would resolve it.

## What this changes

The honest statement about vep-rs is no longer "caching never helps". It is:

> Caching vep-rs loses on transcript-only consequence calling and wins as soon as
> the recipe touches one genome-wide population database. The crossing point is
> measurable, it was predicted from the break-even rule before being measured,
> and it falls between two configurations a practitioner would actually choose
> between.

That is a sharper result than either a uniform win or a uniform loss, and it is
the first point in the matrix where the rule was used predictively rather than
descriptively.

## Not quoting the authors' own comparisons

Unchanged from `../veprs-pilot-grch38-chr22-diagnostic/ANALYSIS.md`: the README
claims 187–343× against Perl VEP, the preprint abstract 78–284×, the headline
F1 of 1.000000000 is *adjusted*, and both authors are Natera employees. None of
it bears on a direct-vs-cached measurement of vep-rs against itself.
