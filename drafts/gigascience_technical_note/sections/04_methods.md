# Methods

## Software and cache design

VCFcache is implemented in Python and delegates normalized BCF processing and
indexing to bcftools/htslib [@danecek2021bcftools]. A blueprint is a
genotype-free set of normalized alleles. Cache construction applies one frozen
annotation command to the blueprint and stores its output with immutable recipe
and parameter snapshots. Sample annotation splits records into cache hits and
misses, applies the same command to misses, and produces a coordinate-sorted,
indexed BCF. Semantic comparison was performed afterward and was never timed.

## Samples and cache strategies

The external WGS benchmark used open KPGP and SGDP GRCh38 genomes and Harvard
PGP GRCh37 genomes. Inputs were restricted to carried, biallelic autosomal SNPs
and short indels and checked for compression, index readability, coordinate
order, contigs, sample identity and completeness. Three genomes per provider
were reserved for cache construction; evaluation comprised 20 KPGP, 20 SGDP
and 12 PGP genomes. PLINK2 KING screening at kinship 0.0884 excluded a
second-degree-or-closer relationship across the construction/evaluation
boundary.

Matched-assay experiments used six KPGP and six SGDP GRCh38 evaluation genomes,
subset to Twist Human Core Exome intervals padded by 125 bp and ACMG SF v3.3
MANE Select coding intervals padded by 20 bp. The principal matched comparison
combined VEP `--everything`, clinical fastVEP and clinical SnpEff measurements
on these identical 12 inputs and the same gnomAD AF ≥1% cache membership,
yielding 108 direct/cached sample-tool-assay pairs. The wider final clinical
campaign used 12 genomes per assembly across Panel, WES and WGS; frozen sample
membership is provided with its reproducer.

Public strategies were assembly-matched VCFcache gnomAD v4.1 joint
exome/genome caches containing sites reaching AF ≥10% or AF ≥1% in any
frequency stratum. Production Zenodo archives were checksum-verified before
use. Provider-specific cohort blueprints contained alleles from the three
reserved genomes; annotations were built independently for every tool and
recipe.
No documented overlap between evaluation cohorts and gnomAD was identified,
although undisclosed sample-level overlap cannot be excluded.

## Annotation recipes and execution

VEP 115.2 ran offline with the assembly-matched Ensembl cache, `--everything`,
`--hgvsg`, buffer size 500,000 and eight forks. fastVEP 0.3.0 used frozen
release binaries, assembly-matched FASTA and transcript caches. Its core recipe
reported consequences and HGVS; clinical and engineering recipes added the
specified supplementary databases. SnpEff 5.4c used `GRCh38.115` or
`GRCh37.75`; the core recipe emitted ANN consequences, whereas the clinical
recipe subsequently used SnpSift to add ClinVar and gnomAD exome AF/AC/AN
fields [@cingolani2012snpeff]. Binary, reference, database and recipe hashes are
archived. These native configurations intentionally represent realistic
within-tool loads; they do not produce identical annotation fields and were not
used for cross-tool performance ranking.

Each unique sample, annotator, recipe, assay and cache condition contributed
one direct and one cached observation; samples, not technical repeats, were the
inferential units. Direct and cached arms used identical inputs, recipes,
databases, storage and CPU allocations and ran in separate Slurm job steps.
End-to-end wall time included normalization, lookup, annotation, output
assembly, compression and indexing, but excluded cache construction and
semantic comparison. Short Panel/WES tasks used an untimed, separate
calibration input before either retained condition to standardize warm state.

## Correctness comparison

The VEP comparator checked variant keys, retained input fields, genotype, the
CSQ definition and all CSQ entries; it canonicalized CSQ and same-locus record
ordering. `HGNC_ID`-only differences were counted but excluded from pass/fail
because VEP 115.2 varies this field with batching [@vep_issue_1959]. The
fastVEP comparator checked every INFO and FORMAT value and relevant header,
canonicalizing only representational ordering. The SnpEff comparator checked
the complete tested annotation level, including ANN, ClinVar and gnomAD fields,
while canonicalizing same-position record and ANN-entry order.

Four SnpEff cells with one warning-only difference were excluded before
summary statistics. Diagnosis used direct SnpEff 5.4c outside all timing: the
target alone, an adjacent insertion followed by the target, and the reversed
pair. No VCFcache process was involved. Exact inputs, outputs and checksums are
archived with the final campaign.

## Controlled annotation-load experiments

Pipeline scaling used held-out GRCh38 sample KPGP-00319 (4,795,706 variants)
and bundled caches with observed hit rates of 80.23% and 90.26%. Six VEP
`--everything` recipes added 0.5, 1, 2, 4, 7 or 10 ms per transcript consequence
through a plugin that emitted no field. Each load had one direct observation
and one cached observation per cache depth.

A separate fastVEP experiment compared core and dense ten-database recipes on
the same held-out WGS under `taskset` limits of 1, 10 or 32 CPUs. The dense
one-CPU exploratory cell was stopped for impractically low throughput and was
not analysed. This controlled cache fixed coverage near 90% to isolate the
interaction between CPU availability and annotation load rather than estimate
biological coverage.

## Cache construction and amortization

Cache construction was excluded from per-sample timing and recorded separately.
For each cohort cache, the measured build duration was divided by paired median
time saved to obtain the number of subsequent genomes required for break-even.
Effective per-sample time after $S$ uses was calculated as
$T_{\mathrm{cached}}+T_{\mathrm{build}}/S$. Public Zenodo caches have no local
annotation-build cost. Cache allele counts, archive and extracted sizes, peak
resident memory and transient disk use were taken from indexed artefacts and
instrumented builds rather than inferred from input size.

## Outcomes and analysis

The primary outcome was speedup
$T_{\mathrm{direct}}/T_{\mathrm{cached}}$; secondary outcomes were relative
runtime, paired time saved, hit rate and cache-build amortization. Serial
throughput was reported descriptively as 86,400 s divided by median wall time.
Cumulative cohort time returned was
$S(T_{\mathrm{direct}}-T_{\mathrm{cached}})-T_{\mathrm{build}}$ and was
converted to compute days. Medians were calculated across genomes. Pooled
confidence intervals used 10,000 sample-level bootstrap resamples stratified
by provider cohort while preserving paired conditions. Modelled and measured
values were labelled separately. Figures were generated in R with ggplot2,
assembled with patchwork, and exported as SVG, PDF and 600-dpi PNG with frozen
source tables and session information.
