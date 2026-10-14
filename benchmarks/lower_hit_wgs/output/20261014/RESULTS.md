# Lower-hit SGDP WGS follow-up, 14 October 2026

The existing measured comparison now includes all 20 SGDP evaluation genomes,
including two Dinka and one Luo genome with 88.82–88.98% coverage. Their VEP
speedups are 6.84–6.96-fold, compared with 4.53-fold for the original Khomani
San genome. The accompanying fastVEP observations use the original core
recipe; they are not the clinical recipe in the matched assay panel.

## Second public Khomani San genome

The newly downloaded `LP6005677-DNA_D03` contains 6,020,947 prepared autosomal
alleles. Of these, 4,723,348 occur in the bundled any-stratum AF >=1% blueprint:
**78.45% coverage and 1,297,599 misses**. The original `LP6005592-DNA_C05`
contains 6,169,350 alleles, 78.98% coverage and 1,296,702 misses.

| Genome | All alleles | QUAL >=30 | SNPs | Indels |
|---|---:|---:|---:|---:|
| Original Khomani San | 78.98% | 79.39% | 78.61% | 80.72% |
| Second Khomani San | 78.45% | 78.93% | 78.09% | 80.13% |
| Dinka SS6004480 | 88.82% | 89.60% | 89.06% | 87.68% |
| Luo LP6005677-DNA_G01 | 88.87% | 89.47% | 88.73% | 89.57% |
| Georgian LP6005442-DNA_B04 | 94.09% | 94.80% | 93.95% | 94.72% |

The lower coverage recurs in both SNPs and indels. QUAL >=30 filtering changes
coverage by less than one percentage point in each genome. The original
Khomani San indel fraction is within the other SGDP genomes' range. These
findings make an isolated problematic sample, excess low-QUAL variants, or an
unusually large indel fraction less convincing explanations for the gap.

All four previously measured whole-autosome hit rates were reproduced exactly
by the new allele-matching screen. The new genome is **coverage-only**: no
annotation speedup, semantic gate, or relatedness clearance has been assigned
to it. It must pass the usual gates before entering performance inference.

## Genotype-confidence control

The bounded control compares all prepared chromosome-20 alleles with alleles
starting within regions confidently genotyped in all five gVCFs (complete GT,
GQ >=20, MIN_DP >=10 if available, otherwise DP >=10; PASS/dot records).
This is a gVCF confidence proxy, not a validated callable-region truth set.
The source GQ/DP/QUAL distributions are also chromosome-20 diagnostics.

The shared regions cover **59,954,567 bases**. The gap persists after this
restriction:

| Genome | All prepared chr20 alleles | Shared confident chr20 regions |
|---|---:|---:|
| Original Khomani San | 69.67% | 71.98% |
| Second Khomani San | 69.30% | 71.60% |
| Dinka SS6004480 | 76.97% | 79.90% |
| Luo LP6005677-DNA_G01 | 76.66% | 79.90% |
| Georgian LP6005442-DNA_B04 | 77.77% | 81.82% |

The Khomani San source records do not have unusually high fractions of low
genotype quality in this control: GQ <20 occurs in 1.94% and 2.20%, versus
2.20–2.97% for the three comparison genomes. DP <10 occurs in 1.55% and 1.79%,
versus 1.74–2.37%. These simple QC measures do not explain the coverage gap.

Absolute chromosome-20 coverage is lower than whole-autosome coverage for all
five genomes. An independent whole-genome streaming intersection, subsequently
restricted to chr20, reproduced the Georgian indexed hit count of 101,257
exactly. The lower regional rates therefore are not an artifact of the indexed
region query in this check. Their cause has not been resolved; chromosome-20
rates must not be substituted for measured whole-genome coverage or used as
genome-wide estimates. The paired chromosome-20 control addresses confidence
filtering on the same chromosome only.

## Interpretation

The replicated low coverage supports a mismatch between the sampled allele
spectrum and the cache's population-frequency selection. This is consistent
with population-specific variation, but does not isolate ancestry as a causal
factor: depth, mapping, source filtering, library preparation and sample
quality can still contribute. No ethnicity-based speedup estimate is justified
from these few genomes. Users should measure their actual cache coverage.

Whole-autosome coverage, quality subsets, source quality and common-region
results are separate tables. The original timing dataset and earlier figure
exports are retained.
