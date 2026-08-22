# Annotation configuration levels

## Why this exists

The manuscript's claim is a break-even rule: caching returns an annotator's
per-variant work but never its start-up cost. Testing it needs points spanning a
range of per-variant cost, and every point must be a configuration somebody
would actually run.

The campaigns as originally run did not satisfy the second condition. VEP ran
its maximal built-in configuration while fastVEP and SnpEff ran their documented
minimums, so part of the measured spread was configuration rather than
annotator. `TODO.md` recorded this as an unresolved, paper-critical gap.

## The levels

**Level 1 — core.** Each annotator's minimal consequence annotation, no
supplementary sources. This is the minimum-work anchor: it is what makes the
start-up-bound regime visible, and it is the only level at which a fast
annotator can be shown to lose time.

**Level 2 — clinical.** Core plus ClinVar and gnomAD allele frequencies. This is
the realistic level.

VEP's offline cache 115 already ships ClinVar 202502, gnomAD exomes v4.1 and
gnomAD genomes v4.1 alongside dbSNP, SIFT, PolyPhen, COSMIC and the regulatory
build. With `--everything` VEP therefore already sits at level 2, and its
existing campaigns are reused unchanged — VEP is not re-run.

| annotator | level 1 | level 2 |
|---|---|---|
| VEP 115.2 | not measured (would need a new cache) | `--everything --hgvsg`, existing |
| fastVEP 0.3.0 | `--hgvs --no-progress`, existing | + ClinVar + gnomAD via `sa-build` / `--sa-dir` |
| SnpEff 5.4c | bare `ANN` | + ClinVar + gnomAD via `SnpSift annotate` |

## Both assemblies

KPGP and SGDP are GRCh38; PGP is GRCh37. PGP carries the weakest
gnomAD-overlap concern of the three cohorts, so it is required in every
comparison rather than optional, and every annotator therefore needs both
assemblies.

The observed hit rates support treating the cohorts as comparable: against the
gnomAD AF >=1% cache the medians are 91.95% (KPGP), 93.04% (SGDP) and 92.66%
(PGP). PGP is a different provider on a different assembly and lands in the same
place, which is what one expects if the rate reflects the fraction of an
individual genome that is common variation rather than any leakage between
cache and evaluation set.

## Cross-tool comparability

`ANNOTATOR_FOLLOWUP_PLAN.md` governs: the question is how much VCFcache reduces
each annotator's runtime, not which annotator is best. Cross-tool wall times are
context, not speedups. The levels therefore do not force identical annotation
content across tools — they place each tool in a configuration a practitioner
would recognise, at two separated levels of per-variant cost.

## Sources

- **ClinVar** — GRCh38 and GRCh37 releases, contigs renamed to UCSC style to
  match every input used here.
- **gnomAD** — v4.1 exome sites for both assemblies; the GRCh37 set is the
  liftover of the same release, so the two assemblies share one gnomAD version.
  Reduced to `AF`, `AC` and `AN`, since the per-population columns multiply size
  without adding lookup work.

## Invocation

Annotators are driven in their documented single-process form. An earlier
SnpEff runner split the input by contig, which gave the direct arm a parallel
speedup the smaller cached arm could not use and so compared two different
execution shapes. Those measurements are retained in the archive as the
optimised-invocation bound, but the documented serial form is primary.
