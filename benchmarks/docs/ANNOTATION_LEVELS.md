# Annotation configuration levels

## Why this exists

The manuscript's claim is a break-even rule: caching returns an annotator's
per-variant work but never its start-up cost. Testing it needs points spanning a
range of per-variant cost, and every point must be a configuration somebody
would actually run.

The configurations intentionally span native workloads rather than pretending
to form a controlled annotator race. VEP, fastVEP and SnpEff cannot produce an
identical minimalist record: their transcript models, consequence and HGVS
implementations, supplementary fields and output schemas differ. The paper
therefore estimates the effect of caching within each tool and treats cross-tool
wall times only as context.

## The levels

**Level 1 — core.** Each annotator's minimal consequence annotation, no
supplementary sources. This is the minimum-work anchor: it is what makes the
start-up-bound regime visible, and it is the only level at which a fast
annotator can be shown to lose time.

**Level 2 — clinical.** Core plus ClinVar and gnomAD allele frequencies. This is
a realistic native configuration for fastVEP and SnpEff, not an assertion of
identical output between tools.

VEP's offline cache 115 ships ClinVar and gnomAD alongside dbSNP, SIFT,
PolyPhen, COSMIC and the regulatory build. `--everything` enables these
efficiently from VEP's prepared cache and remains the real-world
annotation-rich VEP condition. An exploratory narrower recipe
(`--hgvs --check_existing --af_gnomade`) was not promoted to a publication
campaign because it would still not make VEP output or computational work
identical to fastVEP or SnpEff.

| annotator | level 1 | level 2 |
|---|---|---|
| VEP 115.2 | not measured | `--everything --hgvsg`, existing annotation-rich condition |
| fastVEP 0.3.0 | `--hgvs --no-progress`, existing | + ClinVar + gnomAD via `sa-build` / `--sa-dir` |
| SnpEff 5.4c | bare `ANN` | + ClinVar + gnomAD via `SnpSift annotate` |

## Both assemblies

KPGP and SGDP are GRCh38; PGP is GRCh37. PGP carries the weakest
gnomAD-overlap concern of the three cohorts, so it is required in every
comparison rather than optional, and every annotator therefore needs both
assemblies.

The older 52-genome cohort gave median gnomAD AF >=1% hit rates of 91.95%
(KPGP), 93.04% (SGDP) and 92.66% (PGP). The final campaign uses a frozen
12-genome subset per assembly; its PGP selection differs by two genomes and
has a 95.8% median. These are separate estimates, not conflicting results.

## Cross-tool comparability

`ANNOTATOR_FOLLOWUP_PLAN.md` governs: the question is how much VCFcache reduces
each annotator's runtime, not which annotator is best. Every speed-up pairs
direct and cached execution of the same native recipe. Cross-tool wall times are
context only and are not labelled workload-matched or like-for-like. The
annotation-rich VEP and dense fastVEP experiments, the SnpEff core/clinical
contrast and the synthetic-delay sweep jointly test how benefit changes with
pipeline cost.

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
