# Benchmark storage layout

Live filesystem snapshot: 2026-08-03 14:10–14:30 CEST. Values reported by
`df -h` are filesystem values; directory sizes are rounded from `du`. This is
not an OpenStack quota report.

## ITCCcloud_dev: preparation and raw-data authority

`gvbrowse-preproc` (`appuser@10.133.255.21`) is the preparation VM. Its 2-TiB
data volume is mounted at `/mnt/data`: 1.6 TiB used, 250 GiB available (87%).
The 19-GiB root filesystem has 3.3 GiB available.

The authoritative preparation tree is
`/mnt/data/vcfcache_benchmarks` (about 342 GiB):

| Content | Location | Current size |
|---|---|---:|
| Primary/assay raw sources | `sources/` | 34 GiB |
| External KPGP/SGDP/PGP raw sources | `external_wgs/sources/` | 251 GiB |
| Prepared primary/assay sample VCFs | `samples/` | 5.9 GiB |
| Prepared external sample VCFs | `external_wgs/samples/` | 7.1 GiB |
| External normalization work | `external_wgs/work/` | 7.3 GiB |
| Bundled Zenodo VCFcache caches and archives | `bundled_zenodo_caches/` | 17 GiB |
| Bundled Zenodo blueprints and archives | `bundled_zenodo_blueprints/` | 1.3 GiB |
| External three-genome caches | `external_wgs/cohort_caches/` | pending build |

VEP applications and reference caches are under `/mnt/data/apps`. VEP 113 and
115 together occupy about 182 GiB; VEP 115 occupies about 94 GiB and contains
both `115_GRCh37` and `115_GRCh38`. The preparation checkout is
`/mnt/data/home/appuser-projects/vcfcache`.

Raw external downloads remain on this VM. Only normalized samples, QC and
strategy manifests, built cohort caches, and the required bundled caches are
staged to the Slurm cluster.

## ITCCcloud_prod: Slurm execution and durable results

### Head node

`sl-head` has internal address `10.0.0.212` and temporary floating address
`10.133.241.128`. Its preserved 2-TiB data volume is mounted at `/mnt/data`:
484 GiB used and 1.5 TiB available (25% used). Its 19-GiB root filesystem has
15 GiB available.

Important head-node locations:

| Content | Location | Current size |
|---|---|---:|
| Durable Slurm results and task archives | `/mnt/data/slurm-results/` | 222 GiB |
| Staged primary/assay sample VCFs | `/mnt/data/vcfcache_benchmarks/samples/` | 5.9 GiB |
| Bundled VCFcache caches currently staged | `/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/` | 11 GiB |
| Benchmark checkout | `/mnt/data/vcfcache/` | 0.9 GiB |
| VEP 115 cache, GRCh37 plus GRCh38 | `/mnt/data/apps/ensembl-vep/115/` | 94 GiB |

`/mnt/data/slurm-results` is exported over NFSv4. Workers mount the same
physical directory as `/results`; therefore files below worker `/results` do
not consume a second copy on every worker.

The current 222-GiB result tree includes about 100 GiB for the bundled-AF1
primary WGS campaign, 101 GiB for an earlier retained primary campaign, 19 GiB
for an earlier partial campaign, and 2 GiB for the current assay campaign.
These older archives are retained for auditability and have not been deleted.

### Workers

Each disposable worker has a local 500-GiB OpenStack volume mounted at
`/mnt/data`, plus the shared head-node result export at `/results`:

| Worker | Internal IP | Local used | Local available |
|---|---|---:|---:|
| `sl-w1` | `10.0.0.64` | 82 GiB | 385 GiB |
| `sl-w2` | `10.0.0.75` | 69 GiB | 399 GiB |
| `sl-w3` | `10.0.0.168` | 75 GiB | 393 GiB |
| `sl-w4` | `10.0.0.249` | 75 GiB | 393 GiB |
| `sl-w5` | `10.0.0.159` | 94 GiB | 373 GiB |
| `sl-w6` | `10.0.0.167` | 76 GiB | 392 GiB |

Every worker currently has about 5.9 GiB of primary/assay samples, 5.4 GiB of
bundled VCFcache caches, and 25 GiB of VEP 115 GRCh38 data. External staging
adds the 7.1-GiB normalized cohort, GRCh37 VCFcache caches, custom cohort
caches, and the 24-GiB VEP 115 GRCh37 cache. The least-free worker still has
373 GiB available, so this is comfortably within capacity.

Jobs use `/mnt/data/tmp` for worker-local transient run directories. Current
scratch occupancy ranges from about 27 to 53 GiB per worker. Successful task
directories are copied atomically to `/results/campaigns/<campaign-id>/...`;
failed attempts are retained below the campaign's `attempts/` directory.

## Cache types

- VCFcache's ready-made public caches downloaded from Zenodo:
  `/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/` on the preparation VM,
  head, and each worker.
- Cohort-derived three-genome caches:
  `/mnt/data/vcfcache_benchmarks/external_wgs/cohort_caches/`, built on the
  preparation VM and then copied to the same path on the head and workers.
- Controlled hit-rate self-caches:
  `/mnt/data/vcfcache_benchmarks/controlled_runtime/caches/`, prepared on the
  development VM and staged after the external campaign.
- Ensembl VEP's own offline cache, which is distinct from VCFcache:
  `/mnt/data/apps/ensembl-vep/115/cachedir/homo_sapiens/`.
