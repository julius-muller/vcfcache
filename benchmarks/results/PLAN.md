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
| SnpEff core | GRCh38 | 2:55:28 measured | done |
| SnpEff core | GRCh37 | ~3 h | building (job 1000) |
| SnpEff clinical | GRCh38 | ~12.4 h | building (job 995) |
| SnpEff clinical | GRCh37 | ~8-10 h | building (job 999) |
| fastVEP clinical | GRCh38 | ~2 h | annotations built; cache outstanding |
| fastVEP clinical | GRCh37 | ~2 h | annotations building (job 997) |

The partition's `MaxTime` was raised from 12 h to 24 h so the clinical builds fit
serially. That matters: a split build would still produce a valid cache, but it
would be the one artefact in the matrix not made the way the paper says
annotation is run.

## Outstanding

1. Finish the five caches above.
2. Campaigns: 12 GRCh38 genomes and 12 PGP GRCh37 genomes, Panel/WES/WGS,
   direct and cached, at both levels.
3. Fold SnpEff and the clinical level into the manuscript; register the new
   claims in `claims.tsv`.
4. Reconcile the cohort-confounding wording in `sections/03_findings.md` with the
   requirement that PGP appears in every arm.

## Open questions

- `LP6005592-DNA_C05` has 6.17 M variants against 4.6-5.0 M for its siblings and
  a 78.98% hit rate against a 91-93% band. Its indel fraction matches the others,
  so this reads as a genuinely more diverse genome rather than a preparation
  error, and it belongs in the set as the lower bound. Worth a sentence in the
  text rather than a footnote.
- Two GRCh37 clinical throughput runs on the identical slice differed by 14%
  (707.7 s versus 622.6 s). Sizing uses the slower figure; the variance itself
  should be reported rather than averaged away.
