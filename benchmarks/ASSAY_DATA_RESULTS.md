# HPRC, WES, and panel data-preparation results

Preparation completed 2026-07-31 UTC with implementation commit
`f738401ed760`. Large data and machine-readable manifests are under
`/mnt/data/vcfcache_benchmarks`; they are intentionally outside Git.

## Completion gate

`qc/assay_sample_qc.tsv` contains exactly 120 unique cohort/sample rows and all
120 have `PASS` status. Its SHA-256 is
`e16cf66bfd1ede48887f0260574e9bfc485232d4d3f6cccdb9a8f14749dae2b7`.

| Cohort | Samples | Total records | Records/sample, min–median–max | VCF payload |
|---|---:|---:|---:|---:|
| HPRC R2 | 20 | 101,641,659 | 4,740,657–4,859,038–5,851,995 | 1,672,043,545 B |
| Twist core WES | 50 | 1,108,285 | 20,343–21,345–26,292 | 18,449,515 B |
| ACMG SF v3.3 panel | 50 | 12,851 | 196–260–316 | 290,720 B |

Every final file is a sorted, BGZF-compressed, tabix-indexed, single-sample
`.vcf.gz`; no BCF dataset was created. All records are biallelic, nonsymbolic,
and restricted to the intended contigs. The WES outputs have variant records on
chromosome X in all 50 samples. The panel BED includes chromosome X; 41 samples
carry a non-reference panel record there.

## Frozen interval definitions

| Assay | Merged intervals | Target bases | BED SHA-256 |
|---|---:|---:|---|
| Twist Human Core Exome hg38, chr1–22/X | 191,723 | 33,074,111 | `b1c4f29837061526d8524035524cdd745886c767c8f6819a47e4a75ea63c2221` |
| ACMG SF v3.3, Ensembl 115 MANE Select CDS ±20 bp | 1,823 | 423,617 | `11f6395868425252215de90d0ca485e97b3cb2c462d2ec468be009370b571450` |

The panel builder found MANE Select CDS intervals for all 84 frozen ACMG SF
v3.3 genes. The exact source URLs, byte sizes, upstream MD5 values where
available, and local SHA-256 values are in `manifests/assay_sources.tsv`.
Principal source checksums are:

- HPRC v2.0 GRCh38 wave VCF:
  `6316681cb75a418d518c657160e0e700b3252af6e01965d5304d1dccd11b365b`;
- Twist hg38 BED:
  `d7bafeb53f8130b10724359425e900b91ed4d83ae7488c7a7404d96726c86483`;
- Ensembl GRCh38 release 115 GTF:
  `2f8e31578c3aa2f35646927c4a3b3b0dcf0321e57c0ebd3ecc81afcbc836d1a8`;
- 1000 Genomes chrX VCF:
  `eb3db41d7e0e51473d7a88e69e6bb00c189f7606a934c3299294202dec12b569`.

## HPRC selection

The frozen cohort contains 5 AFR, 3 AMR, 4 EAS, 4 EUR, and 4 SAS samples,
spans both sexes in every superpopulation, and has no identity overlap with the
existing 50 1000 Genomes or seven GIAB samples. The exact IDs and metadata are
in `manifests/selected_hprc_r2_samples.tsv`.

HPRC is graph/assembly-derived and therefore remains a separately reported
robustness cohort. It must not be pooled with the read-called WGS samples. The
Twist cohort is an in-silico interval subset, not a wet-lab capture experiment;
its value is the controlled, matched comparison across input sizes.

## Time and storage

The live source download, source hashing, interval construction, two cohort
splits, 100 derived assay files, and all-120-file QC completed in approximately
25 minutes on the current 32-vCPU VM. Streaming BED matching reduced the final
35 WES transformations to about 15 seconds with eight workers.

The complete benchmark root now occupies 48 GiB: sources 34 GiB, work products
5.0 GiB, final samples 5.8 GiB, retained pilots 3.8 GiB, and QC metadata less
than 1 MiB. `/mnt/data` retained 544 GiB free after completion.
