# Plan: the annotation-configuration matrix

Status as of 2026-10-22.

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

VEP's offline cache 115 ships ClinVar, gnomAD, dbSNP, SIFT, PolyPhen, COSMIC and
the regulatory build, so `--everything` is the realistic annotation-rich VEP
condition. VEP, fastVEP and SnpEff do not expose identical transcript models,
consequence/HGVS implementations or output schemas; a narrower VEP command
would standardize some source names but not the actual annotation work. The
matrix therefore spans realistic native configurations and makes inference
within each tool, not from absolute runtime rankings between tools.

## Levels

- **L1 core** — each annotator's minimal consequence annotation. The
  minimum-work anchor that makes the start-up-bound regime visible.
- **L2 clinical** — core plus ClinVar and gnomAD AF.

| tool | L1 core | L2 clinical |
|---|---|---|
| VEP | not run | **done** (`--everything --hgvsg`, both assemblies; annotation-rich) |
| fastVEP | **done** (both assemblies, 52 genomes) | **done** (both assemblies) |
| SnpEff | GRCh38 **done**; GRCh37 **done** | GRCh38 and GRCh37 **done** |
| vep-rs | **done** GRCh38 (jobs 1242, 1251) | **done** via the plugin path (job 1251); not reachable from the converted cache's own data |

`benchmarks/vep_plugins/SyntheticDelay.pm` remains the controlled continuous
sweep of per-variant cost. It models the axis; it does not compare tools.

## Both assemblies, deliberately

PGP is the cohort with the least gnomAD-overlap concern and it is GRCh37, while
KPGP and SGDP are GRCh38. Every tool therefore needs both.

The concern this answers: ~92% of a random genome's variants falling in the
gnomAD AF>=1% blueprint could look tuned. It is not. The older 52-genome cohort
gave 91.95% (KPGP), 93.04% (SGDP) and 92.66% (PGP). The final 12-genome PGP
selection differs by two samples and gives 95.8%; the two estimates are kept
separate and their exact sample lists are frozen.

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
3. ~~Assess a narrow VEP clinical arm.~~ **Closed** — exact three-tool output is
   not natively achievable; exploratory builds and campaigns were canceled
   before timed measurements.
4. Fold the native VEP, fastVEP and SnpEff results into the manuscript without
   presenting cross-tool wall time as a controlled annotator ranking.
5. Reconcile the cohort-confounding wording in `sections/03_findings.md` with the
   requirement that PGP appears in every arm.
6. ~~Decide the vep-rs arm.~~ **Resolved by job 1251, in the opposite direction
   to the pilot's prediction.** The pilot measured transcripts only. With one
   genome-wide gnomAD source vep-rs runs at 31.95 us/variant and crosses even the
   legacy threshold; the clinical level is 40.76 us. Only bare core (7.42 us)
   stays below. A vep-rs campaign therefore measures a predicted **win** at
   clinical level and a predicted **loss** at core — both sides of the rule in one
   annotator, which no other tool in the matrix provides. See
   `raw/archive/veprs-heavy-grch38-chr22-genomewide/ANALYSIS.md`.
7. Decide whether the fused-routing engine enters the manuscript. It halves
   VCFcache's per-variant overhead (11.03 -> 5.55 us) and so moves the threshold
   from ~29.7 to ~15.0 us/variant, which changes the verdict for one
   configuration (`--everything` on transcripts alone, 14.00 us). The figure is
   scaled from a synthetic harness; validating it against a real WGS cache is the
   open step, and `benchmarks/fused_engine/RESULTS.md` names it too. The engine is
   opt-in and legacy remains the default, so every published number stays legacy
   unless this is deliberately changed.

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
returns per-variant work regardless of whether that work is consequence
calculation, database lookup or fixed per-record overhead. The native-tool
results show that behaviour in realistic configurations, while
`SyntheticDelay.pm` supplies the controlled cost axis.

## Open questions

- `hu6664F3` fails concordance on SnpEff at WES and WGS scale (both levels),
  while its panels and every fastVEP assay pass with matching record counts.
  Retained-output job 1205 isolated one warning-only difference at
  `chr2:11273459 C>T`: cache-build SnpEff adds
  `WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS` to two consequences while direct
  WES SnpEff omits it. The four cells remain excluded; no broad exception is
  applied.

- `LP6005592-DNA_C05` has 6.17 M variants against 4.6-5.0 M for its siblings and
  a 78.98% hit rate against a 91-93% band. Its indel fraction matches the others,
  so this reads as a genuinely more diverse genome rather than a preparation
  error, and it belongs in the set as the lower bound. Worth a sentence in the
  text rather than a footnote.
- Two GRCh37 clinical throughput runs on the identical slice differed by 14%
  (707.7 s versus 622.6 s). Sizing uses the slower figure; the variance itself
  should be reported rather than averaged away.

## vep-rs is the minimum-cost end of the axis

The matrix was built to span per-variant cost. vep-rs extends it downward by an
order of magnitude and, in doing so, crosses to the other side of break-even:

| tool and level | per variant | caching at f = 0.9258 |
|---|---|---|
| SnpEff clinical, GRCh38 | ~2370 us | wins by a wide margin |
| SnpEff core, GRCh38 | 558 us | wins |
| vep-rs clinical, genome-wide gnomAD | 40.8 us | wins |
| vep-rs + one genome-wide gnomAD source | 32.0 us | wins, just |
| vep-rs `--everything`, transcripts only | 14.0 us | loses on legacy; **at break-even on fused** |
| vep-rs core | 7.4 us | **loses at every N** |

Break-even needs per-variant cost above 29.7 us, which is VCFcache's own
27.5 us/variant overhead divided by the hit rate. Below that the per-variant
coefficient of `T_direct − T_cached` is negative and the deficit grows with
cohort size, so no amount of data rescues it.

Two things follow for the manuscript. The rule is falsifiable and worth stating
as such: it names a measurable threshold and predicts a direction on both sides
of it. And the honest recommendation for a tool at vep-rs's cost point is not to
cache — which is a better argument for the rule than a table where caching always
wins.

L2 clinical for vep-rs would need either a converter for the Perl cache's
`all_vars.gz`, which the shipped `storable_to_json.pl` does not read, or a full
native build (`vep-cache-builder --include-variations --variation-vcf-dir ...`,
which also wants an indexed genome FASTA we do not hold locally). Neither is
benchmarking work, and the arithmetic above says neither would move vep-rs across
the threshold: adding allele-frequency lookups to a 14 us recipe does not make it
a 30 us recipe.
