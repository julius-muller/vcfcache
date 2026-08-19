# Supplementary Methods

## External WGS selection and quality control

KPGP and SGDP inputs were selected from the DDBJ/NIG public-human-genomes
autosomal reanalysis, and PGP inputs from participant-contributed public VCFs.
Candidate catalogues recorded provider identifier, source URL, assembly,
declared sample identifier and file metadata before download. All preparation
steps were resumable and emitted machine-readable manifests.

Selection was deterministic and performed before timing results were examined.
Three provider-specific genomes were assigned to cache construction. The
evaluation sets contained 20 KPGP, 20 SGDP and 12 PGP genomes. Training and
evaluation membership did not overlap. KPGP and SGDP remained in GRCh38;
PGP remained in GRCh37. No liftover was applied.

Prepared inputs were single-sample, carried-allele VCFs restricted to
biallelic autosomal SNPs and short indels. Quality control tested compression,
index readability, coordinate order, contigs, sample identity, record count
and file integrity. Stable provider identifiers were retained separately from
potentially non-unique VCF sample labels.

Relatedness screening was stratified by assembly. Common SNPs were pruned and
PLINK2 KING robust kinship was evaluated at 0.0884. The design prohibited a
second-degree-or-closer relationship across the cache-training/evaluation
boundary. One related KPGP pair within the evaluation set was retained and
reported because such structure is realistic in an applied cohort.

## Matched assay extension

The assay-by-annotator extension used 12 genomes selected deterministically
from the independent GRCh38 evaluation set before assay timing: six KPGP and
six SGDP. PGP was not included because its GRCh37 configuration would prevent
a matched assembly comparison. Each selected WGS was subset to Twist Human
Core Exome capture intervals padded by 125 bp and to ACMG SF v3.3 MANE Select
coding intervals padded by 20 bp. The resulting WES and Panel inputs therefore
retained the providers and variant-calling pipelines of the independent WGS
cohorts and did not introduce 1000 Genomes-derived performance samples.

Each input was processed once with VEP 115.2 and fastVEP 0.3.0, directly and
through the gnomAD AF ≥1% strategy. The same genomes' frozen WGS observations
completed the matched two-annotator by three-assay design. All cached outputs
were subjected to the same annotator-specific external semantic comparators as
the full WGS campaign. The experiment was designed to identify the sample-size
boundary at which fixed cache-workflow costs are recovered, not to imply that
every assay benefits.

## Public-cache provenance

Only VCFcache caches distributed from production Zenodo were used for the
ready-made public-cache comparison. The GRCh38 VEP 115.2 AF ≥10% and AF ≥1%
archives had DOIs 10.5281/zenodo.18189447 and 10.5281/zenodo.18190046. Their
GRCh37 counterparts had DOIs 10.5281/zenodo.18189051 and
10.5281/zenodo.18189348. Downloads were checked against frozen archive MD5
values before extraction.

The frequency thresholds are any-stratum thresholds: a site entered the
blueprint when at least one reported gnomAD frequency stratum met the cutoff.
The labels therefore do not mean combined-population AF thresholds.
