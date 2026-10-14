# Additional WGS figures: population coverage and benefit

**Figure: Cache coverage and measured annotation savings across the complete
SGDP evaluation set.** (A) All 20 SGDP GRCh38 evaluation genomes, labelled with
source population metadata and sample identifier. Point size indicates the
number of prepared autosomal SNP and short-indel alleles. Coverage uses the
bundled any-stratum gnomAD AF >=1% allele universe. (B) Each genome's measured
speedup relative to direct execution of the same recipe: VEP `--everything`
and core fastVEP. The latter differs from the clinical fastVEP recipe in the
matched assay figure. Dashed curves show `1/(1-f)`, the approximation for
annotation-dominated runtime at hit fraction `f`; they are not fitted curves.
Colours identify the original Khomani San genome, the additional Dinka/Luo
genomes, and the remaining SGDP genomes. Lower-hit genomes were retained in
the original evaluation and are not performance failures. Population labels
do not establish a causal ancestry effect. No new genome has an imputed or
measured annotation speedup in this figure.

**Figure: Allele coverage by variant type and common genotyped regions.**
Chromosome-20 coverage diagnostic for five GRCh38 SGDP genomes, including a second
Khomani San genome. SNP and indel coverage are measured separately. Common
regions are an intersection of chromosome-20 gVCF intervals with complete GT,
GQ >=20 and MIN_DP >=10 when available, otherwise DP >=10. These intervals
are a confidence proxy, not a validated callable-region truth set. Variants
are selected by their start positions. The baseline is all prepared chromosome-20
alleles, so confidence filtering is compared on the same chromosome. This diagnostic is generated only once
both whole-autosome and common-region results are available.

**Figure: Source genotype and variant quality.** Fractions of carried source
chromosome-20 gVCF records with GQ <20, DP <10 or QUAL <30, before publication preparation
decomposes multiallelic records and removes auxiliary fields. Missing values
are separately retained in the source-quality table and do not count as below
threshold. These fractions are diagnostic quality measures, not estimates of
false-positive rates.

**Figure: Whole-autosome coverage by variant type and QUAL threshold.**
All prepared SNP and short-indel alleles compared with the subset at QUAL >=30,
for the same five genomes. This screen covers all autosomes; it is separate
from the bounded chromosome-20 genotype-confidence control.

AF, allele frequency; DP, depth; GQ, genotype quality; QUAL, VCF variant quality;
SGDP, Simons Genome Diversity Project; SNP, single-nucleotide polymorphism;
VCF, variant call format; gVCF, genomic VCF; WGS, whole-genome sequencing.
