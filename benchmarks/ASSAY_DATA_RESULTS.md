# HPRC, WES, and panel data-preparation results

Preparation completed 2026-07-31 UTC and the WES realism correction was made
with implementation commit `fce7aa3`. Large data and machine-readable manifests are under
`/mnt/data/vcfcache_benchmarks`; they are intentionally outside Git.
Publication-facing snapshots, source citations, access conditions, and
checksums are tracked in [`manifests/`](manifests/) and
[SOURCE_PROVENANCE.md](SOURCE_PROVENANCE.md).

## Completion gate

`qc/assay_sample_qc.tsv` contains exactly 170 unique cohort/sample rows and all
170 have `PASS` status. Its SHA-256 is
`c956ec08322ce3b1f4f7bb358fb9f1a5673d155e06acb755892430f23041ced5`.

| Cohort | Samples | Total records | Records/sample, min–median–max | VCF payload |
|---|---:|---:|---:|---:|
| HPRC R2 | 20 | 101,641,659 | 4,740,657–4,859,038–5,851,995 | 1,672,043,545 B |
| Twist capture-like WES, ±125 bp | 50 | 4,010,555 | 73,266–77,056–95,429 | 66,245,371 B |
| Twist strict-target control | 50 | 1,108,285 | 20,343–21,345–26,292 | 18,449,515 B |
| ACMG SF v3.3 panel | 50 | 12,851 | 196–260–316 | 290,720 B |

Every final file is a sorted, BGZF-compressed, tabix-indexed, single-sample
`.vcf.gz`; no BCF dataset was created. All records are biallelic, nonsymbolic,
and restricted to the intended contigs. The WES outputs have variant records on
chromosome X in all 50 samples. The panel BED includes chromosome X; 41 samples
carry a non-reference panel record there.

## Frozen interval definitions

| Assay | Merged intervals | Target bases | BED SHA-256 |
|---|---:|---:|---|
| Twist strict targets, chr1–22/X | 191,723 | 33,074,111 | `b1c4f29837061526d8524035524cdd745886c767c8f6819a47e4a75ea63c2221` |
| Twist capture-like footprint, targets ±125 bp | 165,014 | 78,026,576 | `f007730dc6165a7ca74ea62f4b889d7b4e90be9e52583da56417ffc37f8833e3` |
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

## WES realism calibration

The original strict-target subset had a median of only 21,345 records. It is
preserved under `samples/GRCh38/wes_twist_core_targets` as a useful on-target
control, but it is no longer the primary WES input.

Padding was measured across all 50 matched samples before selecting the new
definition:

| Padding | Footprint | Records/sample, min–median–mean–max |
|---:|---:|---:|
| ±100 bp | 69,657,594 bp | 62,946–66,350–68,947–81,886 |
| ±125 bp | 78,026,576 bp | 73,266–77,056–80,211–95,429 |
| ±150 bp | 86,171,868 bp | 83,396–87,701–91,325–108,790 |

The ±125-bp footprint was selected transparently because its mean is close to
80,000 while its median closely matches the public UCSF NA12878 diagnostic WES
VCF (77,745 sites). A second genuine NA12878 WES VCF from ARUP has 25,066
sites, illustrating substantial laboratory, capture, caller, and filtering
variation. Both CC-BY-licensed GRCh37 files and indexes from Zenodo record
3597727 are checksum-frozen under `sources/wes_calibration/NA12878_GRCh37`.
They are calibration anchors and are not mixed into the GRCh38 timing cohort.

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

The original source download, source hashing, interval construction, two cohort
splits, 100 derived assay files, and all-120-file QC completed in approximately
25 minutes on the current 32-vCPU VM. Streaming BED matching reduced the final
35 WES transformations to about 15 seconds with eight workers.

The realism correction preserved all 50 target-only controls, generated the 50
new capture-like files in under one minute, and expanded the final QC gate to
170 files.

The complete benchmark root now occupies 48 GiB: sources 34 GiB, work products
5.0 GiB, final samples 5.8 GiB, retained pilots 3.8 GiB, and QC metadata less
than 1 MiB. `/mnt/data` retained 544 GiB free after completion.
