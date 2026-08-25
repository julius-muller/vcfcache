# Final clinical campaign reproduction record

This directory freezes the code and inputs that produced the accepted
2026-08-24 campaign. Files ending in `.frozen.*` are byte-for-byte copies from
the Slurm head after jobs 1102, 1105, and 1163--1165 completed. Their checksums
are in `SHA256SUMS`.

Important limitations of the frozen runner are documented rather than hidden:

- parallel array tasks split each logical archive across a base and `-job<id>`
  directory;
- only array task zero wrote `RUN_MANIFEST.json`, so GRCh38 manifests recorded
  KPGP but omitted SGDP;
- GRCh37 samples were discovered with `ls | head -12`.

The canonical non-destructive correction is
`../../results/raw/final_campaign/RUN_MANIFEST.json`, and `samples.tsv` freezes
the exact cohort. Future campaigns must use the explicit sample manifest and a
single controller-created output directory.

`diagnose_hu6664f3.sh` is an untimed retained-output diagnostic. It is not part
of the publication timing campaign.

`diagnose_snpeff_context.sh` and `diagnose_snpeff_minimal.sh` are untimed,
direct-SnpEff-only follow-ups. Jobs 1236 and 1238 proved that the one SnpEff
warning difference is processing-order dependent and reproduced it with two
records without invoking VCFcache. Results and the minimal inputs are archived
under `../../results/raw/diagnostics/hu6664F3/snpeff_context/`.

## Aborted narrow VEP exploration

The `vep_clinical_*`, `build_vep_clinical.sh`, `smoke_vep_clinical.sh`,
`vep_compare_cli.py` and `vep_matched_campaign.sh` files preserve an exploratory
2026-08-25 design. Its smoke tests passed, but cache-build jobs 1226 and 1231 and
dependent campaign jobs 1232 and 1233 were canceled before any timed campaign
task started. These files did not produce the accepted 2026-08-24 campaign and
must not be used as evidence or figure input. The design was abandoned because
VEP, fastVEP and SnpEff cannot natively produce identical transcript
consequences, HGVS representations or output schemas; matching a few source
names would not create a controlled cross-tool workload.
