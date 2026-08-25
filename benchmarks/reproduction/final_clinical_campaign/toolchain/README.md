# Frozen toolchain

Captured read-only on worker `sl-w2` by Slurm job 1214. `versions.txt` records
the runtime versions and VCFcache commit. `assets.sha256` and
`assets.bytes.tsv` freeze the fastVEP binary, FASTA/index and transcript caches;
all supplementary fastVEP `.osa2` databases; SnpEff/SnpSift binaries,
configuration and predictor databases; the ClinVar and gnomAD VCFs/indexes;
bcftools; and the VEP container.

The VEP transcript cache is a directory, so `vep_cache_files.sha256` contains a
sorted checksum for every file and `vep_cache_manifest.sha256` freezes that
complete manifest.

The campaign ran VCFcache 0.5.3 at commit
`cd871861c3b18b25315df560fcaf374b2decb527`; the later manuscript checkout may
have a newer version and must not be substituted in provenance.
