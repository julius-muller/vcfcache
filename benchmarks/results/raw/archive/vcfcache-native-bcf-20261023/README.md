# Resident native BCF pilot, 2026-10-23

New engineering experiment on `feature/native-bcf-routing`, based on
`8830bf7b3`. Existing publication measurements and figures are unchanged.

Primary experiment:

- `fastvep-1273`: two real SGDP WES inputs, complete genome-wide cache first,
  then a Twist-region standard-BCF projection. Fused/native/direct, one run each.
- `veprs-core-1276`: the same experiment with vep-rs core. All audits passed.
- `veprs-clinical-1277`: stopped after direct disagreed with the two cached
  engines at 22/82,714 records. **Not a validated clinical speedup.**
- `clinical-diagnostic-1281`: two discrepant loci rerun in isolation; those
  outputs agree with the full direct run, not the cache-derived annotations.
- `final-1285`/`final-1286`: final resident streaming path checked against saved
  direct outputs on the same workers; no new direct baseline was rerun.
- `controls-1287`: three alternating paired repetitions on an explicitly
  perfect-hit control retaining real WES locus distribution, then a WGS
  buffer-fallback smoke check. Not a biological hit-rate finding.
- `calibration-1288`/`1289`/`1290`: direct CLI runs on misses and one variant per
  occupied contig, for fastVEP, vep-rs core and vep-rs clinical respectively.
  Input preparation is outside timing and on local worker storage.

All conditions use separate Slurm steps, 16 allocated CPUs, local staged cache
BCFs, and complete-record/relevant-header audits outside timing. Transcript and
supplementary annotation resources retain their frozen existing paths. This is
an engineering pilot, not a randomized publication campaign. Cross-node timings
and distinct frozen source snapshots must not be pooled.

Cache provenance:

- fastVEP: existing locally built clinical cache derived from the bundled
  GRCh38 gnomAD AF>=1% blueprint, not a bundled fastVEP download.
- vep-rs: newly built core/clinical BCF caches from the same bundled blueprint.
  Build provenance retains original durations and freezes resource identities.
  Full GRCh38 transcripts were converted from existing VEP cache resources:
  job 1269, 3,102/3,102 files, 552,279 transcripts, no conversion warnings.
- `measured_sources`: frozen A/B/C deployment sources. A is an invalid pilot
  (bad indexed access); B is the primary experiment; C streams final recipe
  merging while releasing retained records. Final repository sources additionally
  fix recipe-returned-record counters for incomplete annotations and conservatively
  replay case-insensitive allele matches. Executable hashes distinguish builds
  sharing the development 0.3.0 label.
- Setup scripts are archived separately, including the original cache builder
  and explicit staged-only metadata repair. The source cache recipes and true
  build duration were not overwritten.

Failed preparation/pilot attempts remain visible in logs: 1271 empty-directory
setup, 1270 output-index defect, 1274 invalid recipe metadata. Jobs 1275,
1282–1283 were cancelled dependencies. Calibration/control setup 1278–1280 and
1284 was retried with indexed intermediate inputs; these failed before timed
conditions. No failed condition is quietly turned into a successful comparison.

Large BCFs, transcript JSON, cache payloads and scratch data remain on the
cluster. They are not committed here. See the maintained
[results report](../../../../native_bcf/RESULTS.md) for interpretation.
