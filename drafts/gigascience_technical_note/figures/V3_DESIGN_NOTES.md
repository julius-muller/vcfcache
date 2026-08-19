# V3 figure design and evidence boundary

The v3 assets use new filenames and do not overwrite the v2 figure set.

## Figure 1

Figure 1 combines one-time setup and per-sample execution in one connected
schematic. Panel A distinguishes three supported starting points:

1. a ready-to-use, recipe-specific cache downloaded from Zenodo;
2. a Zenodo blueprint annotated once with the user's fixed recipe; or
3. variants from an existing cohort, converted to a blueprint and annotated
   once.

The recipe-specific cache feeds panel B directly. Semantic colour roles are
consistent: blue denotes downloaded/lookup assets, orange denotes user inputs,
purple denotes annotation work, green denotes reusable annotations and final
cache/output products, and grey denotes transformations or compatibility
checks. Branch labels use title case and describe actions rather than slogans.

## Figure 2

The historical Panel/WES campaign is not eligible for the v3 main figure
because its samples were derived from 1000 Genomes and therefore overlap the
source universe of gnomAD. A matched extension instead uses six KPGP and six
SGDP GRCh38 genomes from the independent external evaluation set. Twist Human
Core Exome capture intervals and ACMG SF v3.3 MANE Select coding intervals are
applied to the same genomes. VEP and fastVEP are measured once per sample and
assay against the AF ≥1% cache strategy; their already-frozen WGS rows complete
the matched 2 annotators x 3 assays x 12 genomes design.

An initial extension attempt was excluded after short-assay timings showed a
large cold-state order effect. In one WES example, cached-first/direct-second
times were 181.6/27.3 seconds, whereas a direct-first/cached-second example was
261.5/116.8 seconds. The corrected campaign runs one separate KPGP-00319 assay
subset through the cached path before either measured condition in every task.
That warm-up is untimed, excluded from all summaries and distinct from the 12
evaluation genomes. Exactly one direct and one cached measurement remain for
each sample/tool/assay cell.

The corrected extension completed all 48 new Panel/WES cells without an
execution or semantic failure. Joined to the 24 matched frozen WGS cells, it
produced 72 validated observations. Median Panel/WES/WGS speedups were
0.71/2.29/9.50 for VEP and 0.31/0.54/2.02 for fastVEP. Figure 2 therefore
shows the measured break-even boundary rather than implying that caching is
advantageous for every assay and recipe.

PGP is intentionally absent from this matched main panel because its GRCh37
configuration changes the direct fastVEP baseline. In the complete external
WGS data, median direct fastVEP runtime was 207 seconds for PGP versus 402--427
seconds for the two GRCh38 cohorts, while fixed VCFcache work did not halve.
The resulting lower PGP relative speedup is therefore an assembly/recipe and
fixed-overhead interaction, not evidence of poorer biological cache coverage.
All 52 genomes remain visible in Supplementary Figure 1.

## Figure 3

Panel A replaces the confusing direct-annotation identity series with arrows
from direct to cached runtime. Arrow length is the observed wall time returned
per genome. Panel B retains relative speedup but adds the theoretical ceilings
calculated from each observed hit rate:

$$
\lim_{T_{\mathrm{direct}} \to \infty} \mathrm{speedup}
= \frac{1}{1-f}.
$$

Thus absolute hours returned continue to increase with pipeline cost while the
relative speedup approaches a finite hit-rate ceiling. The two observations
are complementary rather than contradictory.
