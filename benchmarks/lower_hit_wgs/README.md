# Lower-hit real-WGS follow-up

This additive analysis uses the original 20 SGDP evaluation genomes and a
coverage-only screen of one additional public Khomani San genome. It does not
alter the matched 12-genome Figure 2 dataset or replace earlier figures.

## Existing measured comparison

The RStudio workspace at
`/home/j380r/pCloudDrive/workspace/jumu_rscripts/vcfcache` contains a new
`LOWER-HIT WGS` section in `figures_gigascience.R`. It uses only the existing
52-genome campaign's VEP `--everything` and **core** fastVEP measurements.
These fastVEP rows must not be described as the clinical recipe in Figure 2.

The figure includes all 20 SGDP genomes rather than selecting only low-hit
examples. All use GRCh38 and the bundled any-stratum gnomAD AF >=1% allele
universe. Panel A shows hit rate and variant burden; Panel B shows measured
within-recipe speedup and the annotation-dominated approximation `1/(1-f)`.
Population labels are source metadata, not independently inferred ancestry.

Existing low-coverage genomes are:

| Sample | Population label | Hit rate | VEP speedup |
|---|---|---:|---:|
| LP6005592-DNA_C05 | Khomani San | 78.98% | 4.53x |
| SS6004480 | Dinka | 88.82% | 6.96x |
| LP6005677-DNA_G01 | Luo | 88.87% | 6.84x |
| LP6005443-DNA_H08 | Dinka | 88.98% | 6.85x |

The archived preparation QC gives an indel fraction of 17.62% for the original
Khomani San genome, within the 16.8–17.9% range of these 20 SGDP genomes. Its
lower total hit rate therefore cannot be explained simply by a larger indel
fraction. SNP/indel-specific hit rates and source call quality require the
variant-level screen below.

## New untimed screen

`screen_wgs.py` downloads `LP6005677-DNA_D03` from the same public DDBJ/NIG
SGDP reanalysis, reuses the publication's `prepare_one()` normalization
implementation, and preserves all new artifacts in an isolated directory.
Hom-reference and missing-genotype gVCF records are excluded before that helper
because they cannot contribute a carried allele. A chr1:1,000,000–3,000,000
control using the original Khomani San source confirmed identical sets of
5,567 normalized allele/genotype keys with and without this early exclusion
(SHA-256 `16ab631f4fc148dc22a3cbf0b48fab63374bb22f72458475b970a812dac381cf`).
The complete source remains available for source-quality and region analysis.
It compares the two Khomani San genomes, one Dinka, one Luo, and a Georgian
genome with high cache coverage. No new annotator timings are generated.

The analyses are:

- Whole-autosome allele coverage, separately for SNPs and indels.
- The same analysis restricted to QUAL >=30 variants.
- Source gVCF carried-record GQ/DP/QUAL diagnostic counts on chromosome 20
  before FORMAT fields are removed by benchmark preparation.
- A bounded chromosome-20 intersection of high-confidence genotyped regions
  across all five sources: PASS/dot records, complete GT, GQ >=20, MIN_DP >=10 when
  present or DP >=10 otherwise. The coverage comparison uses variant-start
  membership in these zero-based half-open intervals. It is compared with all
  prepared chromosome-20 alleles, not with whole-autosome coverage, to avoid
  confounding chromosome choice and confidence filtering.

The common-region mask is a **gVCF-derived proxy**, not an independently
validated callable-region truth set. It does not establish equal false-positive
rates or isolate a causal ancestry effect. Source QC counts are gVCF records;
coverage counts are normalized carried biallelic alleles. Whole-autosome
SNP/indel and QUAL screens are genome-wide; genotype-confidence controls are
explicitly restricted to chromosome 20 and cannot establish genome-wide QC.

The selected cache blueprint came from the bundled Zenodo artifact. Its exact
allele matching is measured with `bcftools isec -c none`. This estimates cache
membership, not semantic equivalence of newly annotated outputs. A new genome
must pass relatedness and semantic gates before entering performance inference.

Example invocation on the preparation VM:

```sh
python benchmarks/lower_hit_wgs/screen_wgs.py \
  --root /mnt/data/vcfcache_benchmarks/external_wgs \
  --output /mnt/data/vcfcache_benchmarks/lower_hit_wgs_20261014 \
  --reference /mnt/data/resources/reference/ucsc/hg38.fa.gz \
  --cache /mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
```

Source metadata and technical background:

- [Primary SGDP sample table](https://discovery.ucl.ac.uk/id/eprint/1516771/12/Balloux_5_18_2016_SGDP_supplementary_data_table_1_sample_information_table.xlsx).
- [Second Khomani San public source](https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/SGDP/GVCF/LP6005677-DNA_D03.hs38DH/).
- [SGDP manuscript](https://doi.org/10.1038/nature18964).
- [DDBJ reanalysis methods](https://doi.org/10.1093/nar/gkaf1273).

Results and captions are recorded in the dated output directory after collection.
