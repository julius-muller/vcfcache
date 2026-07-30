# Public benchmark VCF setup

## Output

Large data live under `/mnt/data/vcfcache_benchmarks`:

```text
sources/{1000g,giab}/       immutable downloaded files
samples/GRCh38/{1000g,giab} final per-sample VCFs and TBI indexes
samples/GRCh38/all          flat symlink view of all 57 VCF/index pairs
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
