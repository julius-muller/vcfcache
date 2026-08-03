# Publication benchmark implementation plan

## Headline claim

VCFcache changes annotation cost from scaling with all input variants to scaling
mainly with variants absent from the cache, while preserving annotation-equivalent
output.

## Implementation phases

1. **Preparation test suite** — implemented
   - Test deterministic sample selection, source-manifest parsing, resumability,
     and VCF validation before downloading the full cohort.
   - Gate full preparation on local tests, a live-network smoke test, and disk
     preflight.
2. **Public sample preparation** — implemented
   - Prepare 50 autosomal 1000 Genomes VCFs and seven GIAB VCFs as described in
     [DATA_SETUP.md](DATA_SETUP.md).
   - Freeze sample IDs, source checksums, transformations, and final checksums.
2a. **Assay-size and representation-robustness cohorts** — implemented
   - Test the HPRC selection, ACMG gene definition, MANE BED construction,
     chrX completion, interval subsetting, final format, and QC before bulk
     preparation.
   - Add 20 HPRC Release 2 graph-derived GRCh38 samples as a separately
     interpreted robustness cohort.
   - Derive 50 matched Twist-core WES and 50 matched ACMG SF v3.3 panel VCFs
     from the primary identities, including chromosome X without mutating the
     established autosomal WGS files.
   - Keep every final dataset as sorted, BGZF-compressed, tabix-indexed VCF;
     freeze source and output checksums.
   - Retain strict Twist targets as a mechanism control and use a separately
     labeled ±125-bp, empirically calibrated capture-like WES cohort for the
     primary size comparison.
   - Completed with 20 HPRC, 50 capture-like WES, 50 target-only WES, and 50
     matched panel files. See [ASSAY_DATA_RESULTS.md](ASSAY_DATA_RESULTS.md).
3. **Benchmark runner** — implemented; pilot complete
   - Pin the VCFcache commit, container, bcftools, VEP 115.2 cache, CPUs, RAM,
     scratch storage, and thread counts.
   - Capture cgroup CPU, memory peak, I/O, wall time, variant counts, hit rate,
     cache size, and output concordance.
4. **Headline experiment** — paired pilot complete; Slurm cohort in progress
   - The confirmed HG02079 pilot achieved 9.56× end-to-end speedup at a 96.01%
     cache-hit rate with semantic equivalence across all 4,329,621 variants.
     See [PILOT_RESULTS.md](PILOT_RESULTS.md).
   - The memory-heavy final statistics scan was replaced with streaming
     counting and the exact pair was repeated at commit `88c018086b21`.
   - Run cached and uncached VEP 115.2 `--everything` on 49 primary WGS
     samples, 20 capture-like WES samples, and 20 ACMG-panel samples.
   - Use one validated paired measurement per sample, with condition order
     balanced across samples. Samples are the independent units; a full-cohort
     warm-up and technical triplicates do not add independent evidence.
   - Retain a small repeatability control only: the observed HG02079 duplicate
     differed by 0.84% uncached and 0.08% cached, far below between-sample
     variation.
5. **Mechanism and complexity experiments** — implemented; run pending
   - Use deterministic self-caches on one real capture-like WES input to measure
     observed hit rates near 50%, 80%, 90%, 95%, and 100%.
   - Compare vanilla VEP, VEP `--everything`, and vanilla VEP with deterministic
     no-output delays of 5 or 20 ms per transcript consequence.
   - Run four uncached baselines and 20 cached cells once each. This technical
     mechanism experiment requires no repeats; every cell must pass the same
     semantic comparison as the real-cohort runs.
   - Use the 49 WGS, 20 WES, and 20 panel pairs for the real-input-size curve;
     do not spend the publication window on redundant short-run repeats.
6. **Analysis and manuscript** — implemented; final data pending
   - Treat samples as the statistical unit and report paired ratios, medians,
     IQRs, and bootstrap 95% confidence intervals.
   - Publish raw metrics, manifests, environment locks, plotting code, model
     JSON, and figure source data.
   - Generate separate assay, independent-WGS/cache-strategy, controlled-model,
     and simplified user-impact SVGs directly from the archived tables.

## Main performance figure

- **A — Avoided work:** representative WGS flow and stacked timing breakdown,
  showing cache hits, misses, and variants sent to VEP.
- **B — Real WGS performance:** paired cached/uncached runtimes across 50 samples,
  colored by superpopulation, with speedup and CPU-hours saved.
- **C — Mechanism:** measured speedup versus controlled cache-hit rate, faceted
  by annotation complexity and overlaid with the fitted runtime model.
- **D — Cohort economics:** cumulative CPU-hours for 1–10,000 genomes, showing
  public-cache and custom-cache scenarios plus the observed break-even point.

## Infrastructure

Data preparation and sequential pilots run on the current 32-vCPU, 62-GiB VM.
The 50-sample repeated annotation matrix should now move to Slurm. After the
streaming fix, GNU time's largest-process peak RSS is 5.22 GiB uncached and
749 MiB cached; observed aggregate VEP worker RSS was approximately 20–23 GiB
during the uncached run.
