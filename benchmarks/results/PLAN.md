# Plan: the annotation-configuration matrix

Status as of 2026-08-22.

## Why this matrix exists

The manuscript's central claim is a break-even rule: caching returns an
annotator's per-variant work but never its start-up cost, so whether it helps
depends on how much per-variant work the recipe performs. Testing that needs
points spanning a range of per-variant cost, and each point has to be a
configuration somebody would actually run.

The earlier matrix failed the second test. VEP was benchmarked realistically
(`--everything --hgvsg`), but fastVEP and SnpEff were benchmarked with **zero
supplementary databases**. The headline spread was therefore partly an artefact
of configuration.

## The fact that shapes it

VEP's offline cache 115 already ships ClinVar 202502, gnomAD exomes and genomes
v4.1, dbSNP 156, SIFT, PolyPhen, COSMIC and the regulatory build. With
`--everything`, **VEP is already at the clinical level** — it needs no rerun.

That fixes the target for the other two: adding ClinVar + gnomAD allele
frequencies brings fastVEP and SnpEff to a level VEP already occupies, using
sources all three support natively.

## Levels

- **L1 core** — each annotator's minimal consequence annotation. The
  minimum-work anchor that makes the start-up-bound regime visible.
- **L2 clinical** — core plus ClinVar and gnomAD AF.

| tool | L1 core | L2 clinical |
|---|---|---|
| VEP | not run — would need a new VEP cache (~10 h) | **done** (`--everything`, both assemblies) |
| fastVEP | **done** (both assemblies, 52 genomes) | annotations built; caches outstanding |
| SnpEff | GRCh38 **done**; GRCh37 building | GRCh38 and GRCh37 building |

`benchmarks/vep_plugins/SyntheticDelay.pm` remains the controlled continuous
sweep of per-variant cost. It models the axis; it does not compare tools.

## Both assemblies, deliberately

PGP is the cohort with the least gnomAD-overlap concern and it is GRCh37, while
KPGP and SGDP are GRCh38. Every tool therefore needs both.

The concern this answers: ~92% of a random genome's variants falling in the
gnomAD AF>=1% blueprint could look tuned. It is not. Hit rates are 91.95% (KPGP),
93.04% (SGDP) and 92.66% (PGP) — PGP is a different provider on a different
assembly and lands in the same place, which argues the rate is biological rather
than a property of the sample selection.

## Cache builds

Immutable per recipe, so each level needs its own genome-wide build over the
18.87 M-allele blueprint.

| cache | assembly | estimate | state |
|---|---|---|---|
| SnpEff core | GRCh38 | 2:55:28 | **done** |
| SnpEff core | GRCh37 | 3:01:38 | **done** |
| SnpEff clinical | GRCh38 | 13:11:20 | **done** |
| SnpEff clinical | GRCh37 | 7:20:53 | **done** |
| fastVEP clinical | GRCh38 | 28:49 | **done** |
| fastVEP clinical | GRCh37 | 6:05 | **done** |

The two fastVEP estimates were ~2 h; both were wrong by more than an order of
magnitude, in the direction of fastVEP being far cheaper than scaling from
SnpEff suggested. Estimates in this table are labelled as such.

The partition's `MaxTime` was raised from 12 h to 24 h so the clinical builds fit
serially. That was necessary rather than precautionary: the GRCh38 clinical build
took **13:11:20** and would have been killed at the old limit. That matters: a split build would still produce a valid cache, but it
would be the one artefact in the matrix not made the way the paper says
annotation is run.

## Outstanding

1. ~~Build the five caches.~~ **All five are built and content-verified.**
2. ~~Campaigns.~~ **Done** — 180 cells, 176 concordant. See `RESULTS.md`.
3. Fold SnpEff and the clinical level into the manuscript; register the new
   claims in `claims.tsv`. **Next step.**
4. Reconcile the cohort-confounding wording in `sections/03_findings.md` with the
   requirement that PGP appears in every arm.

## Assembly is a cost axis for one annotator and not the other

Measured on the built caches, chr22, annotations per record:

| | GRCh37 | GRCh38 | ratio | build time ratio |
|---|---|---|---|---|
| fastVEP | 3.19 | 12.14 | 3.8x | 4.7x (6:05 -> 28:49) |
| SnpEff | 4.27 | 12.29 | 2.9x | 0.97x (3:01:38 -> 2:55:28) |

Both annotators agree the assemblies differ: `GRCh37.75` is Ensembl release 75
and carries a far smaller transcript set than `GRCh38.115`, so every variant has
roughly a third of the consequences to compute. fastVEP's runtime tracks that
almost exactly. **SnpEff's does not move at all** -- it took marginally longer on
the assembly with a third of the work.

So SnpEff's per-variant cost is dominated by fixed per-record overhead rather
than per-transcript work, while fastVEP's is dominated by the annotation itself.

Two consequences for how the results are read:

- **Hit rates stay comparable across assemblies.** Cache membership is
  annotator-independent, so the 91.95 / 93.04 / 92.66% across KPGP / SGDP / PGP
  remains a fair comparison and still defends the ~92% figure.
- **Speedups are not comparable across assemblies for fastVEP.** A lower fastVEP
  speedup on the GRCh37 PGP genomes would mean "less per-variant work to return",
  not a cohort effect. Presenting them side by side without that caveat invites
  exactly the wrong reading. For SnpEff the assembly is close to a null variable
  and PGP can be read as a straight replicate.

This strengthens the break-even argument rather than complicating it: caching
returns per-variant work regardless of whether that work is real computation or
fixed overhead, and here are two tools at the same nominal annotation level that
differ in which one it is. That is a less contrived demonstration than the
`SyntheticDelay.pm` sweep, because nothing was arranged to produce it.

## Open questions

- `hu6664F3` fails concordance on SnpEff (all assays, both levels) while passing
  on fastVEP, with matching record counts. Unexplained; suspected batch-dependent
  non-determinism in SnpEff's `ANN`. One genome of twelve.

- `LP6005592-DNA_C05` has 6.17 M variants against 4.6-5.0 M for its siblings and
  a 78.98% hit rate against a 91-93% band. Its indel fraction matches the others,
  so this reads as a genuinely more diverse genome rather than a preparation
  error, and it belongs in the set as the lower bound. Worth a sentence in the
  text rather than a footnote.
- Two GRCh37 clinical throughput runs on the identical slice differed by 14%
  (707.7 s versus 622.6 s). Sizing uses the slower figure; the variance itself
  should be reported rather than averaged away.
