# Publication source provenance

This document is the human-readable provenance record for the VCFcache
publication benchmark. The machine-readable companion is
[`manifests/source_catalog.tsv`](manifests/source_catalog.tsv). The source
inventory and recorded checksums were audited on 2026-07-31 UTC.

## Provenance hierarchy

The study uses four layers of records:

1. This document identifies each resource, its scientific role, release,
   coordinate system, citation, and data-use status.
2. Tracked files in [`manifests/`](manifests/) freeze URLs, upstream and local
   checksums, selected identities, and derived interval definitions.
3. Operational manifests under `/mnt/data/vcfcache_benchmarks` contain the
   absolute VM paths and SHA-256 of every final VCF.
4. The eventual publication archive must contain the two QC tables, raw timing
   tables, environment records, plotting inputs, and the final Git commit.

Checksums identify bytes, not merely filenames. Upstream MD5 values validate
downloads when providers publish them; SHA-256 values identify local and
derived artifacts. The local UCSC FASTA is BGZF-recompressed, so its compressed
checksum intentionally differs from the vendor gzip while the decompressed
streams are identical (SHA-256
`056974f68c0e2d006cd76b3c1218d0e6ede326653e7d92518cc7ad72ff11920c`).

## Source inventory

| Resource | Frozen release and coordinate system | Role | Detailed machine record |
|---|---|---|---|
| 1000 Genomes | NYGC 30× high-coverage phased panel, 2022; GRCh38 | Primary read-called WGS cohort; chr1–22 and assay-only chrX | `source_files.tsv`, `assay_sources.tsv` |
| Population metadata | Phase 3 integrated call sample panel, 2013 | Population, superpopulation, and provider `gender` field used for deterministic balancing | `selected_1000g_samples.tsv` |
| Genome in a Bottle | NIST GIAB v4.2.1; GRCh38; HG001–HG007 | External benchmark genomes | `source_files.tsv` |
| HPRC | Release 2/v2.0 minigraph-cactus `wave.vcf.gz`; GRCh38 coordinates | Graph/assembly-derived robustness cohort | `assay_sources.tsv`, `selected_hprc_r2_samples.tsv` |
| KPGP | DDBJ/NIG public-human-genomes autosomal gVCFs; GRCh38/hs38DH | Korean real-world WGS training and held-out evaluation | generated `external_wgs_candidates.tsv`, selection, and QC manifests |
| SGDP | DDBJ/NIG public-human-genomes autosomal gVCFs plus Mallick et al. sample table; GRCh38/hs38DH | Region-stratified real-world WGS after documented HGDP/1000G exclusion | generated external-WGS manifests and KING report |
| Harvard PGP | Participant-contributed public WGS VCFs; native GRCh37/hg19 | Provider-diverse real-world WGS; downloaded header must establish GRCh37 and exactly one sample | generated external-WGS manifests; no phenotype metadata retained |
| Twist | Human Core Exome covered-target BED for hg38, file dated 2024-09 | Strict WES target control and padded capture-like footprint | `assay_sources.tsv`, `assay_regions.tsv` |
| Ensembl | Release 115 human GTF; GRCh38 | MANE Select CDS coordinates | `assay_sources.tsv`, `assay_regions.tsv` |
| ACMG | Secondary Findings v3.3 Table 1 | Membership of the 84-gene small panel | `config/acmg_sf_v3.3_genes.txt` |
| Real WES anchors | NA12878 ARUP and UCSF VCFs, Zenodo 3597727; GRCh37 | Variant-count calibration only; not performance inputs | `assay_sources.tsv` |
| Reference sequence | UCSC GRCh38.p14 soft-masked FASTA, GCA_000001405.29 | Contig bounds when padding BED intervals | `source_catalog.tsv` |
| Public cache universe | gnomAD v4.1 joint exome+genome sites; native GRCh38 plus published GRCh37 liftover | Assembly-matched variants eligible for pre-annotation | `source_catalog.tsv` and cache artifact hashes |

The complete direct URLs are deliberately stored in TSV rather than abbreviated
in manuscript prose. The 1000 Genomes table records all 22 autosomal VCF/index
pairs and all seven GIAB VCF/index pairs. The assay table records HPRC, Twist,
Ensembl, chrX, and both real-WES VCF/index pairs. Thus a changed `latest`
symlink or vendor page cannot silently change the frozen study identity.

## Public-cache derivation

The benchmark cache is not an opaque third-party artifact. It was derived from
the public gnomAD v4.1 joint exome+genome sites Hail Table at
`gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht`
using [`scripts/export_gnomad_hail.py`](../scripts/export_gnomad_hail.py).
The resulting VCF header records Hail 0.2.137-733ac4ccd943.

Variants were retained when **any** element of `joint.freq` had AF ≥ 0.01.
The exported INFO/AF and INFO/AC values came from the combined `joint.freq[0]`
element. The resulting sites-only BCF contained 18,869,857 records. VCFcache
removed genotypes, split multiallelic records, and built an unannotated
blueprint with the same record count. Ensembl VEP 115.2 then annotated that
blueprint using offline cache 115, GRCh38, `--hgvsg`, `--everything`, buffer
500,000, eight forks, and eight bcftools threads.

Consequently, artifact labels must say **any-stratum AF ≥ 1%**. A site can be
included because one ancestry, sex, or subset stratum passes 1% even when its
combined AF is below 1%. This distinction affects expected WGS coverage and is
not recoverable from the exported INFO/AF field alone.

gnomAD v3.1 incorporated approximately 2,500 1000 Genomes samples into its
HGDP/1KG subset, and the v3.1.2 genomes were carried into v4. The study's NYGC
WGS samples and the selected HPRC identities are therefore not an independent
cohort for estimating this cache's real-world WGS hit rate. They remain valid
for runtime calibration conditional on observed hit rate and for disjoint
cohort-derived cache experiments.

## Bundled-cache provenance gate

The publication public-cache matrix contains only cache bundles discoverable
through VCFcache and downloaded from production Zenodo. The allow-list and
archive checksums are frozen in `manifests/bundled_caches.tsv`:

| Strategy | Shipped alias | Zenodo DOI | Archive MD5 |
|---|---|---|---|
| any-stratum AF ≥ 10% | `cache-gnomad-v4.1-GRCh38-joint-af01-vep115.2-e` | `10.5281/zenodo.18189447` | `088cf426461a51b77bdfcd5dcd2233f4` |
| any-stratum AF ≥ 1% | `cache-gnomad-v4.1-GRCh38-joint-af001-vep115.2-e` | `10.5281/zenodo.18190046` | `3ac438461eac0cf42c75717156d7b2d4` |
| any-stratum AF ≥ 10%, GRCh37 | `cache-gnomad-v4.1-GRCh37-joint-af01-vep115.2-e` | `10.5281/zenodo.18189051` | `96bb1edd0e067d9c933256bd112e4589` |
| any-stratum AF ≥ 1%, GRCh37 | `cache-gnomad-v4.1-GRCh37-joint-af001-vep115.2-e` | `10.5281/zenodo.18189348` | `f7d246a7adf44b778d6dc1383153eff2` |

The runner retains the downloaded archives and writes a validated provenance
record beside each extracted cache. Missing or conflicting provenance stops
the benchmark. Zenodo record `18190666` is a blueprint rather than an annotated
cache bundle and is excluded from this ready-made-cache arm. Cohort-derived
caches are separately labelled custom and never represented as bundled. The
external runner resolves these allow-listed bundles by sample assembly and
rejects cross-assembly cache use.

## Official-blueprint provenance gate

Alternative annotation scenarios use the separately verified official GRCh38
blueprints, not an undocumented reconstruction. Their frozen input identities
are in `manifests/bundled_blueprints.tsv`:

| Strategy | Shipped blueprint alias | Zenodo DOI | Archive MD5 |
|---|---|---|---|
| any-stratum AF ≥ 10% | `bp-gnomad-v4.1-GRCh38-joint-af01` | `10.5281/zenodo.18190424` | `c3d1ea67acd62b3fd9f30ea132d98a41` |
| any-stratum AF ≥ 1% | `bp-gnomad-v4.1-GRCh38-joint-af001` | `10.5281/zenodo.18190436` | `6b7403ff03815500ba49c52ad285746c` |
| any-stratum AF ≥ 0.1% | `bp-gnomad-v4.1-GRCh38-joint-af0001` | `10.5281/zenodo.18190666` | `1e44e7c08c8fb6aec6913eb2914ffabc` |

The runner retains each archive, verifies both byte size and Zenodo MD5, and
writes a role-specific provenance record beside the extracted blueprint. A
cache built from one of these inputs is reproducible and publication-eligible
as a local annotation scenario, but is not itself a bundled cache. Its
annotation YAML, parameter snapshot, output cache, and build duration become
additional derived-artifact provenance.

This gate was exercised on 2026-08-02 on the preparation VM and independently
on the Slurm head. Identical extracted blueprint files and provenance records
are available on `sl-w1` through `sl-w6`; the Zenodo archives are deliberately
excluded from worker staging. The cluster was idle and its queue empty after
staging.

Artifact identities are:

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| Derived gnomAD AF≥0.01 BCF | 257,844,944 | `46c4441095c09f91bc38832d6afba36172c3d2d3e963a7001c37c4c56412a727` |
| VCFcache blueprint | 102,913,821 | `2828116ca9a996fbd5259ea7cb477ec4499a2c01c77c6662d2c345610695057b` |
| VEP-annotated cache | 3,950,678,136 | `e9ae4a550c42a7cb1e68ccc1ef358c5e3558dab998c7bf833faa4f0571e1da7e` |
| VEP Apptainer SIF | 238,407,680 | `3adb1a1c49a43c425aadc3b9fc00c162144625cafe4a6759fe0a5e8d564fb114` |
| Frozen `annotation.yaml` | 931 | `cda583eab6436aea3b18ec6ac7768a9062a9d115aa4bedc38f5f3e7dda385973` |
| Frozen `params.snapshot.yaml` | 813 | `02a14bb11b80074a549b5611c19a27d1a90e15dc77e5b90441b2ff03c9a40f8a` |

The first comment in the immutable annotation snapshot says “GRCh37” because
the file was adapted from that recipe. This is a stale comment only: its
executable `genome_build`, the substituted VEP `-a` argument, cache requirement
check, contig dictionary, and VEP cache metadata all specify GRCh38. The
snapshot must remain unchanged; this clarification should accompany the
archive.

## Reference-genome identity

The interval builder used `/mnt/data/resources/reference/ucsc/hg38.fa.gz`, an
indexed BGZF copy of the UCSC GRCh38.p14 soft-masked FASTA. Its 711 sequences
sum to 3,299,210,039 bp, matching UCSC's p14 chromosome-size file. UCSC reports
MD5 `efb4edf237dd3594e94610ed803c8a44` for the provider gzip. The local BGZF is
SHA-256 `68dabbf6c07654f5ffcda4975afb48ceb233ca549b3459dcc5c7e7289d9fa4b9`;
both yield the identical decompressed SHA-256 reported above. Only contig
lengths were read from the `.fai`; bases were not used to call or normalize
study variants.

## Data use, attribution, and redistribution

- 1000 Genomes/IGSR data are public. Use the project terminology guidance when
  describing populations; this study treats labels as source metadata and
  reports performance strata, not biological racial categories.
- HPRC states that its data are in the public domain. The HPRC Data Use Best
  Practices still require appropriate attribution and responsible population
  descriptions.
- NIST encourages publication use of GIAB data without embargo and requests
  citation of the relevant dataset README and benchmark papers.
- Ensembl project-generated data are available without restriction, while
  third-party content can retain separate terms.
- The NA12878 WES benchmark record is licensed CC BY 4.0.
- Twist makes the target BED available for download and use. Because the vendor
  page does not assign a standalone open-data license, cite the page and
  checksum, and recheck its terms before redistributing the BED in a paper
  deposit. The derived interval checksum can always be published.
- The ACMG gene symbols are a reproducibility record derived from a cited policy
  table; the policy article itself is not redistributed.

Only openly accessible, de-identified public resources are used. No new human
participants were recruited and no attempt is made to re-identify individuals.
The manuscript should use the institution's standard determination language
for whether analysis of public de-identified data requires ethics review.

## Citations to retain

- Byrska-Bishop *et al.* High-coverage whole-genome sequencing of the expanded
  1000 Genomes Project cohort including 602 trios. *Cell* (2022).
  https://doi.org/10.1016/j.cell.2022.08.004
- Wagner *et al.* Benchmarking challenging small variants with linked and long
  reads. *Cell Genomics* (2022). https://doi.org/10.1016/j.xgen.2022.100128
- Lucas *et al.* HPRC2: A human pangenome reference with near-complete coverage
  of common genetic variation. *bioRxiv* (2026).
  https://doi.org/10.64898/2026.07.21.739710
- MANE Collaboration. A joint NCBI and EMBL-EBI transcript set for clinical
  genomics and research. *Nature* (2022).
  https://doi.org/10.1038/s41586-022-04558-8
- Dyer *et al.* Ensembl 2025. *Nucleic Acids Research* (2025).
  https://doi.org/10.1093/nar/gkae1071
- ACMG Secondary Findings Working Group. ACMG SF v3.3 list. *Genetics in
  Medicine* (2025). https://doi.org/10.1016/j.gim.2025.101454
- NA12878 WES Benchmark dataset. Zenodo (2019).
  https://doi.org/10.5281/zenodo.3597727
- gnomAD Production Team. gnomAD v4.1 release and joint frequency resource
  (2024). https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1
- gnomAD Production Team. gnomAD v3.1 HGDP and 1000 Genomes subset
  (2020). https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/
- Genome Reference Consortium. Human reference assembly GRCh38.p14,
  GCA_000001405.29. UCSC download and checksums:
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p14/

Provider landing pages and current data-use statements are linked from
[`MATERIALS_AND_METHODS.md`](MATERIALS_AND_METHODS.md), because those webpages
should be cited with an access date in addition to the archival papers.

## Submission-time provenance gate

Before manuscript submission, all of the following must be true:

- the final benchmark Git commit and VCFcache version are inserted into Methods;
- the complete full-cohort `environment.json` and scheduler/cgroup records are
  archived, not only the single-sample pilot;
- `qc/sample_qc.tsv`, `qc/assay_sample_qc.tsv`, raw run metrics, and figure source
  data are deposited with a DOI and linked from this repository;
- the VEP SIF recipe or immutable image identifier and the VEP cache install log
  are archived alongside their checksums;
- the final Slurm allocation, CPU model, RAM request, storage type, OS, replicate
  randomization seed, and excluded/failed runs are reported; and
- Twist redistribution terms are rechecked; if redistribution is not clearly
  allowed, publish the URL, source checksum, transformation code, and derived
  checksum rather than the vendor BED itself.
