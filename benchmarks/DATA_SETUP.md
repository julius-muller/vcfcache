# Public benchmark VCF setup

## Output

Large data live under `/mnt/data/vcfcache_benchmarks`:

```text
sources/{1000g,giab}/       immutable primary-cohort downloads
sources/{hprc_r2,twist,ensembl}/ assay and robustness sources
samples/GRCh38/{1000g,giab} primary per-sample VCFs and TBI indexes
samples/GRCh38/hprc_r2      20 HPRC R2 robustness VCF/index pairs
samples/GRCh38/wes_twist_core 50 matched WES VCF/index pairs
samples/GRCh38/panel_acmg_sf_v3.3 50 matched panel VCF/index pairs
samples/GRCh38/all          flat symlink view of all 57 VCF/index pairs
regions/GRCh38/             frozen merged WES and panel BED files
manifests/                  selected samples, URLs, checksums, and provenance
qc/                         per-sample statistics and validation reports
work/                       resumable partial downloads and chromosome shards
logs/                       command logs
```

The root filesystem is not used. `TMPDIR` is set below the benchmark root.

## 1000 Genomes

Source: NYGC 30x high-coverage, integrated phased GRCh38 panel (2022).

- Select 10 unrelated samples from each of AFR, AMR, EAS, EUR, and SAS.
- Select five male and five female samples per superpopulation.
- Rank eligible samples by SHA-256 of
  `vcfcache-paper-v1:<sample_id>` for deterministic selection.
- Exclude the seven GIAB identities before selection.
- Use autosomes 1–22.
- Keep SNPs and indels carried by the selected individual.
- Preserve cohort-level `AF`, `AC`, and `AN` without recalculation.
- Retain `GT` as the only FORMAT field.
- Split all selected samples in one source pass per chromosome.
- Concatenate chromosomes, coordinate-sort, BGZF-compress, and tabix-index.

## Genome in a Bottle

Source: NIST GIAB GRCh38 v4.2.1 small-variant benchmark VCFs for HG001–HG007.
Use the archived HG002 v4.2.1 file so all seven samples use the same release.
Retain the official individual-level INFO and FORMAT fields, then sort,
BGZF-compress, and index the final copy.

## Validation

Every final file must:

- pass `bgzip -t`;
- have a valid `.tbi`;
- contain exactly one sample;
- contain records only on chromosomes 1–22;
- be coordinate sorted;
- have a recorded variant count and SHA-256 checksum.

1000 Genomes files additionally must contain only SNPs/indels, only non-reference
sample genotypes, INFO tags `AF`, `AC`, and `AN`, and FORMAT tag `GT`.

The final `qc/sample_qc.tsv` records cohort, population, sex, counts, file size,
checksum, and validation status for all 57 VCFs.

## HPRC Release 2 robustness cohort

Source: the official HPRC v2.0 minigraph-cactus wave VCF in GRCh38
coordinates. The source is graph/assembly-derived and is analyzed as a
separate robustness cohort, not pooled statistically with read-called 1000
Genomes samples.

- Intersect the 232 VCF sample columns with official 1000 Genomes Phase 3
  population metadata.
- Exclude all identities already represented by the 50-sample primary cohort
  or seven GIAB genomes.
- Deterministically choose 20 samples with seed
  `vcfcache-paper-hprc-r2-v1`: AFR 5, AMR 3, EAS 4, EUR 4, and SAS 4.
  There are only three eligible AMR samples with the frozen metadata, so all
  three are retained and the two remaining slots go to AFR.
- Use chromosomes 1–22, split multiallelic sites, keep carried biallelic SNPs
  and indels with REF and ALT shorter than 50 bp, and retain GT plus source
  AF/AC/AN.
- Split all 20 samples in one source pass, then emit ordinary VCF.gz/TBI pairs.

## Matched WES cohort

Source intervals: Twist Human Core Exome covered targets for hg38, downloaded
from Twist's current data-file page. Freeze the vendor file checksum, restrict
to chromosomes 1–22 and X, sort, and merge overlapping/bookended intervals.

Subset the same 50 primary 1000 Genomes identities. Autosomal variants come
from the immutable prepared WGS inputs; chromosome X comes from the matching
official NYGC 30x chrX callset. The output is an *in-silico captured WES VCF*,
which controls ancestry, genotype, and caller while changing only target size.
It must not be described as a wet-lab exome callset.

## Matched small-panel cohort

Gene definition: the 84 symbols in ACMG SF v3.3 Table 1. The exact list is
frozen in `benchmarks/config/acmg_sf_v3.3_genes.txt`. Generate GRCh38 targets
from Ensembl release 115 GTF records satisfying all of:

- feature is `CDS`;
- `gene_name` is in the frozen 84-gene list;
- transcript has the `MANE_Select` tag; and
- contig is chromosome 1–22 or X.

Convert GTF coordinates to BED, pad each CDS by 20 bp on both sides, then sort
and merge. Assert that every one of the 84 genes contributes at least one
interval. Subset the same 50 samples and chrX shards used for WES. This cohort
is intended to expose the cache/startup crossover at small input sizes, even
if cached execution is not faster.

## Assay-cohort validation

`qc/assay_sample_qc.tsv` must contain exactly 20 HPRC, 50 WES, and 50 panel
rows. Every final VCF must pass BGZF testing, have a valid tabix index, contain
exactly one sample, use only the allowed contigs, contain no symbolic or
multiallelic records, and have nonzero indexed record count, byte size, and
SHA-256 recorded. Source URLs, upstream MD5 where available, local SHA-256,
BED definitions, interval counts, and target bases are frozen under
`manifests/`.
