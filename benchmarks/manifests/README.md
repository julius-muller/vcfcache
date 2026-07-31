# Frozen publication manifests

These small, tracked files are the manuscript-facing snapshot of the much
larger machine manifests under `/mnt/data/vcfcache_benchmarks/manifests`.
They make the cohort definition and source identity reviewable without access
to the preparation VM.

- `source_catalog.tsv` is the resource-level provenance index, including
  release, assembly, role, source URI, checksums, access date, and reuse note.
- `source_files.tsv` freezes every 1000 Genomes autosomal and GIAB VCF/index
  URL and the upstream MD5 when one was published.
- `assay_sources.tsv` freezes downloaded HPRC, Twist, Ensembl, chrX, and real
  WES calibration artifacts. Paths are relative to the benchmark root.
- `selected_1000g_samples.tsv` and `selected_hprc_r2_samples.tsv` freeze all
  selected identities and the deterministic selection seeds.
- `assay_regions.tsv` freezes the three derived BED definitions and their
  SHA-256 checksums. Paths are relative to the benchmark root.

The large operational manifests remain authoritative for individual output
artifacts: `qc/sample_qc.tsv` contains every WGS/GIAB output checksum and
`qc/assay_sample_qc.tsv` contains all 170 HPRC/WES/panel output checksums.
Before archival, copy those two QC tables and the raw timing tables into the
paper's data deposit and record the deposit DOI here.

`accessed_utc` records when a source was downloaded or reverified for this
study. An empty checksum means that the provider did not publish that checksum;
it does not mean that the local artifact was left unhashed.

