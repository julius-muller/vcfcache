# Run log

What was actually run on the cluster, in order. One line per job. Tool output
itself stays on the cluster (`/mnt/data/slurm-results`); this is the index to it.

Cluster: 6-node `benchmark` partition, `sl-w1`–`sl-w6`, 16 vCPU / 30 GB each,
`OverSubscribe=EXCLUSIVE`. Jobs run as `appuser`. No Slurm accounting, which is
why this log exists.

## 2026-08-21 — SnpEff core, GRCh38

| job | what | outcome |
|---|---|---|
| 945 | SnpEff core Panel+WES, 12 genomes, **split** invocation | superseded — split gives the direct arm a parallel speedup the cached arm cannot use |
| 966 | SnpEff core WGS, 3 genomes, split | superseded; retained as the optimised-invocation bound |
| 974 | SnpEff core WGS, 12 genomes, **serial** | primary — median 8.17x, median hit 92.58% |

`MaxTime` raised 12 h → 24 h on the partition (live + `slurm.conf`); needed, as
the GRCh38 clinical build later ran 13:11:20.

## 2026-08-22 — gates and cache builds

| job | what | outcome |
|---|---|---|
| 983/984 | clinical gate, first attempt | failed — `SnpSift annotate` needs `-c snpEff.config` |
| 990/991 | clinical gate, second | failed — a `grep -qm1` guard tripped `pipefail` via SIGPIPE |
| 992 | clinical gate GRCh38 | **PASSED**; 327 records/s |
| 993/994 | clinical gate GRCh37 | 994 **PASSED** after fixing `genome_build` in params; 422 records/s |
| 996 | core gate GRCh37 | **PASSED** — validates the newly installed `GRCh37.75` |
| 997/998 | fastVEP ClinVar + gnomAD `.osa2`, both assemblies | built |
| 995 | cache `snpeff_clinical_grch38_af001` | 13:11:20 |
| 999 | cache `snpeff_clinical_grch37_af001` | 7:20:53 |
| 1000 | cache `snpeff_core_grch37_af001` | 3:01:38 |
| 1001/1002 | caches `fastvep_clinical_grch38/37_af001` | gated, then 28:49 / 6:05 |

All five caches content-verified: full allele counts, ANN/CSQ + ClinVar + gnomAD
present, cross-checked between annotators on chr22.

## 2026-08-23 — campaigns v1, invalidated

| job | what | outcome |
|---|---|---|
| 1003, 1010–1013 | first campaigns | **excluded** — stale `c.bcf` across assay iterations with an unchecked return code produced fabricated 0.5 s runtimes |
| 1021–1023 | GRCh37 Panel/WES, after lifting the intervals with UCSC `liftOver` | **excluded** — same defect |
| 1099–1101, 1152 | diagnostics | identified all three defects |
| 1102–1106 | campaigns v2 | GRCh38 arms primary; GRCh37 still failed on a contig-order digest |
| 1163–1165 | campaigns v3 | GRCh37 arms primary |

Three defects, all in the harness, none in VCFcache: stale outputs, duplicate
input records in five PGP genomes, and a digest assuming a shared contig order.
The user-facing consequences of the latter two are in `/KNOWN_LIMITATIONS.md`.

## 2026-08-24 — final

180 cells complete, 176 concordant. Results in `RESULTS.md`, all rows in
`campaign_metrics.tsv`.

Primary runs, by archive directory:

- `*-grch38-v2-serial*` — fastVEP clinical, SnpEff clinical (GRCh38)
- `*-grch37-v3-serial*` — fastVEP clinical, SnpEff core, SnpEff clinical (GRCh37)

## 2026-08-25 — audit corrections and aborted narrow VEP exploration

| job | what | outcome |
|---|---|---|
| 1199 | retained-output `hu6664F3` diagnostic | failed before comparison — submitted as `ubuntu`, which cannot execute the `appuser` VCFcache environment |
| 1201 | diagnostic under `appuser` | annotation completed; comparator failed because core SnpEff correctly has no clinical INFO tags |
| 1205 | corrected core + clinical WES diagnostic | **complete** — one warning-only difference at `chr2:11273459 C>T` in both levels; exact values archived under `raw/diagnostics/hu6664F3/` |
| 1214 | toolchain fingerprint | **complete** — binaries, versions, references, transcript caches and supplementary databases frozen under `benchmarks/reproduction/final_clinical_campaign/toolchain/` |
| 1203 | first matched-VEP recipe smoke, both assemblies | passed, then superseded because it included gnomAD genome and exome frequencies while the fastVEP/SnpEff level uses the prepared exome source |
| 1206 | first VEP cache-build attempt | failed before annotation — build script pre-created an invalid empty database directory |
| 1210 | second VEP cache-build attempt | failed before annotation — VEP statistics directory absent |
| 1215 | third VEP cache-build attempt | failed before annotation — VEP HTML statistics path not container-visible; annotator reports removed consistently |
| 1220 | broader-load cache build | canceled when the exome-only load mismatch was identified; no measurements accepted |
| 1225 | exact ClinVar + gnomAD-exome VEP smoke, both assemblies | **passed** |
| 1226 | exact-load builds | GRCh37 task completed only by detecting/reusing the wrong broader cache and was excluded; GRCh38 task canceled at 09:03 CEST before completion; no measurements accepted |
| 1231 | exact-load GRCh37 rebuild with recipe-identity guard | canceled by design at 09:03 CEST before completion; no measurements accepted |
| 1232/1233 | dependent narrow-load VEP campaigns, GRCh37/GRCh38 | canceled while still dependency-pending; no array task started |
| 1236 | direct-SnpEff context matrix for `chr2:11273459 C>T` | **complete** — target alone is warning-free from both sample and blueprint; the blueprint ±1 kb context introduces the warning |
| 1237 | first minimal direct-SnpEff reproducer | failed before annotation because of a shell-local initialization error; no result produced |
| 1238 | corrected minimal direct-SnpEff reproducer | **complete** — preceding `chr2:11273458 A>AC` adds the warning; reversing the same two records removes it |
| 1241 | compact public-ticket reproducer verification | **complete** — all three minimal VCFs reproduced the expected warning pattern and target ANN digests |

The narrow-load campaign was abandoned after review established that common
source labels would not produce identical annotation work or output across VEP,
fastVEP and SnpEff. Existing VEP `--everything` measurements remain the native
annotation-rich VEP result. The exploratory scripts are retained as provenance
only and are excluded from summaries, figures and claims. The recipe-identity
guard remains a useful general safeguard for future cache builds.

## 2026-10-22 — vep-rs 0.4.0, pilot

| job | what | outcome |
|---|---|---|
| 1242 | vep-rs pilot: convert Perl VEP 115 chr22 cache to JSON, then four configurations at two input sizes | **complete** — diagnostic; see `raw/archive/veprs-pilot-grch38-chr22-diagnostic/` |
| — | field-capability probe on the converted cache, head node (`V2_veprs_field_probe.sh`) | **complete** — no clinical fields reachable |

The converted cache carries transcripts only, so the reachable ceiling is core
plus HGVS; `--everything` adds no SIFT, PolyPhen, gnomAD AF or ClinVar, and
warns that `--regulatory` has no effect. Two-point fit: 5.68 µs/variant (core)
to 14.00 µs/variant (`--everything`), start-up 0.16–0.25 s. Against VCFcache's
21.7 s + 27.5 µs/variant, break-even needs ≥ 29.7 µs/variant at *f* = 0.9258, so
the rule predicts caching **loses at every scale** for this tool. That negative
prediction, not a speedup, is what a vep-rs arm would test.

## 2026-10-22 — vep-rs annotation-heavy regime, and the fused-routing engine

| job | what | outcome |
|---|---|---|
| 1243, 1245, 1247, 1249 | heavy-ladder attempts | failed before measuring; causes in `V5_heavy_ladder.sh` header |
| 1250 | heavy ladder, partial | `L0`/`L1` measured, then died on a `printf '%s'` that left `\t` literal in a TSV header |
| 1251 | heavy ladder, genome-wide gnomAD | see `raw/archive/veprs-heavy-*/` |
| — | fused-routing threshold check (local) | `raw/archive/fused-routing-threshold-check/` |

**One genome-wide gnomAD source already crosses break-even.** Job 1250 measured
vep-rs core at 7.49 µs/variant and core plus `gnomADc` over the 99.9 M-record
genome-wide joint gnomAD at **34.4 µs/variant**, with the field asserted present
on 3,540,863 of 3,787,554 output rows. Against the 29.7 µs legacy threshold that
is a win on a single supplementary source — not the two-to-three the
exome-restricted probe had suggested, because the exome file simply has no AF to
return off-target.

**Fused routing roughly halves the threshold.** Per-variant overhead 11.03 µs
legacy against 5.55 µs fused (1.99×), while start-up moves only 0.752 → 0.729 s.
That is the right signature: caching returns per-variant work and never start-up.
Rescaling the production 27.5 µs/variant gives **~15.0 µs/variant**, down from
29.7. vep-rs `--everything` at 14.0 µs then sits essentially *at* break-even on
transcripts alone. Estimate scaled from the branch's synthetic harness, not a
production measurement; validating it on a real WGS cache is the open step.
