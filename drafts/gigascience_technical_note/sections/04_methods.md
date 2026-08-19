# Methods

## Software design

VCFcache is implemented in Python and delegates variant transformation,
compressed BCF I/O and indexing to bcftools/htslib [@danecek2021bcftools]. A
blueprint is a genotype-free collection of normalized variant alleles. Cache
construction applies one annotation command recipe to that blueprint and
stores the annotated records together with immutable annotation and parameter
snapshots. At sample annotation time, records are split into hits and misses
against an indexed cache, misses pass through the configured annotation
command, and records are merged into a coordinate-sorted indexed BCF.

The command recipe receives explicit input and output paths and exposes the
same integration boundary for annotators that consume VCF streams. This paper
benchmarks VEP and fastVEP. Other command recipes are supported by the public
interface but are not assigned performance claims here.

Semantic comparison was an external post-processing step and did not
contribute to reported wall time.

## Real-world WGS datasets

The external benchmark used open autosomal WGS resources from KPGP, SGDP and
Harvard PGP. KPGP and SGDP gVCFs were obtained from the DDBJ/NIG public human
genomes reanalysis and retained in GRCh38/hs38DH coordinates. SGDP metadata
were used to exclude documented HGDP/1000 Genomes identities and balance
geographic regions [@mallick2016sgdp]. Participant-contributed PGP VCFs were
accepted only when their headers established GRCh37 coordinates and exactly
one sample [@ball2012pgp]. Numeric PGP contigs were renamed to UCSC form
without coordinate liftover. No phenotype fields were collected.

Selection, download, preparation and quality control were deterministic and
recorded in repository manifests. Prepared inputs retained carried, biallelic
autosomal SNPs and short indels. Each file was required to pass BGZF integrity,
tabix-index readability, coordinate-order, contig, single-sample and non-zero
content checks. Broad record-count bounds were used only to detect truncated
or misclassified inputs, not to select samples by a performance outcome.

Three genomes per provider cohort were assigned to cache construction and
were excluded from evaluation. Evaluation comprised 20 KPGP, 20 SGDP and 12
PGP genomes. PLINK2 KING screening was stratified by assembly. Kinship above
0.0884 across the cache-training/evaluation boundary failed preparation. One
second-degree KPGP pair was retained within the evaluation set to preserve
realistic cohort structure.

For the matched assay extension, six KPGP and six SGDP GRCh38 evaluation
genomes were selected deterministically before assay timing. Each was subset
to Twist Human Core Exome capture intervals padded by 125 bp and ACMG SF v3.3
MANE Select coding intervals padded by 20 bp. Their frozen WGS measurements
were retained, yielding matched Panel, WES and WGS inputs without introducing
1000 Genomes-derived performance samples.

## Cache strategies and annotation recipes

The public strategies were the VCFcache-distributed, assembly-matched gnomAD
v4.1 joint exome/genome caches for sites reaching AF ≥10% or AF ≥1% in any
frequency stratum. GRCh38 annotated cache records were downloaded from Zenodo
records 18189447 and 18190046; the corresponding GRCh37 records were from
18189051 and 18189348. Archive checksums were verified before extraction and
execution. Three-genome provider-specific blueprints supplied the cohort-cache
strategy. Blueprint membership was shared between annotators, but its
annotations were rebuilt separately with the assembly-matched VEP or fastVEP
recipe.

VEP 115.2 ran offline with its assembly-matched Ensembl cache, --everything,
--hgvsg, a buffer size of 500,000 and eight forks. VEP output was converted to
BCF with eight bcftools threads. fastVEP 0.3.0 used frozen release-mode
binaries, assembly-matched reference FASTA and transcript caches. Exact
binary, reference and recipe fingerprints are retained in the campaign
manifests. The manuscript archive will identify VCFcache [FINAL RELEASE],
commit [FINAL BENCHMARK COMMIT], its container digest and all external asset
checksums.

## Benchmark execution

Each external-WGS sample was one Slurm task with one common direct baseline and
three cached conditions. The matched assay extension used one direct condition
and the AF ≥1% cached condition for each sample, assay and annotator. Conditions
ran as separate Slurm job steps, providing
condition-specific process accounting and peak-memory measurement. Execution
order was balanced deterministically across samples. Every unique sample,
annotator and cache condition was timed once; samples, not technical repeats,
were the inferential units. A small repeatability control was kept separate
from the inferential dataset.

Short Panel/WES cells received one untimed pre-run of a separate KPGP-00319
assay subset before either retained condition. This standardized warm
annotator, transcript and cache-lookup state and was excluded from reported
wall time. The calibration input was not among the 12 assay-extension genomes
and no retained sample/tool/assay cell was repeated.

The cached and direct arms used identical input, annotator, recipe, databases,
storage, CPU allocation and output contract. End-to-end wall time included
input processing, cache lookup, annotation, output assembly, compression,
and indexing. Cache construction and semantic comparison were excluded.
Resource records and logs were archived atomically per task, and failed
attempts were retained rather than silently replaced.

## Correctness comparison

The VEP comparator traversed cached and direct indexed BCFs in the same
canonical contig order. It compared CHROM, POS, REF and ALT; retained input
AF, AC and AN values; genotype; the CSQ definition; and the complete set of
CSQ entries. It canonicalized CSQ order and complete-record order within
identical CHROM:POS loci. Differences confined to HGNC_ID were counted and
reported but excluded from pass/fail because VEP 115.2 can vary this field
with buffering and batching, as documented in Ensembl VEP issue 1959
[@vep_issue_1959]. Any
other difference failed the pair.

The fastVEP comparator tested every INFO and FORMAT value and relevant header
definition. It canonicalized only INFO-tag order, CSQ-entry order and
complete-record order within identical loci. It ignored no annotation field,
including FV_CLINVAR, FV_GNOMAD, FV_1KG and FV_TOPMED when present. Comparator
reports were generated outside timed cells.

## Controlled pipeline-complexity experiment

Pipeline scaling used held-out sample KPGP-00319, a GRCh38 WGS with 4,795,706
variants. The AF ≥10% and AF ≥1% bundled caches produced 80.23% and 90.26%
hits, respectively. Six VEP --everything recipes invoked a custom plugin that
slept for 0.5, 1, 2, 4, 7 or 10 ms per transcript consequence and returned no
annotation fields. The delay applied only to variants submitted to VEP.

Each load had one direct observation and one cached observation per cache
depth, with alternating order and no technical repeats. The AF ≥1% extension
reused the frozen direct baselines from the same sample and recipes. All
cached outputs were required to match the designated direct reference.

## Outcomes and statistical analysis

The primary outcome was relative wall time
$R_{\mathrm{cached}} = T_{\mathrm{cached}}/T_{\mathrm{direct}}$. Secondary
outcomes were speedup $1/R_{\mathrm{cached}}$, paired wall time saved, cache
hit rate and absolute time returned. Medians were estimated across genomes.
Pooled 95% bootstrap intervals used 10,000 sample-level resamples stratified by
KPGP, SGDP and PGP, preserving all paired conditions for each resampled genome.
Deterministic seeds are recorded in the R plotting scripts.

The prespecified runtime model was
$T_{\mathrm{cached}} \approx T_{\mathrm{overhead}} +
(1-f)T_{\mathrm{direct}}$. For the controlled complexity
experiment, fixed overhead was estimated from the measured zero-delay
configuration, and predicted relative runtime was compared with each delayed
measurement. Cache-build amortization was calculated as
$T_{\mathrm{effective}} \approx T_{\mathrm{cached}} +
T_{\mathrm{build}}/S$, where $S$ is the number of cache uses. Modelled values
are labelled separately from direct measurements.

## Figure generation

All plots were generated in R with ggplot2. Composite figures were assembled
with patchwork and exported as SVG, PDF and 600-dpi PNG. Frozen source tables,
R session information and renderer records accompany every final figure
snapshot.
