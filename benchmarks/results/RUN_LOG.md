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
