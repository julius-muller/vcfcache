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
2a. **Assay-size and independent robustness cohorts** — implemented
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
4. **Headline experiment** — paired pilot complete; cohort pending
   - The confirmed HG02079 pilot achieved 9.56× end-to-end speedup at a 96.01%
     cache-hit rate with semantic equivalence across all 4,329,621 variants.
     See [PILOT_RESULTS.md](PILOT_RESULTS.md).
   - The memory-heavy final statistics scan was replaced with streaming
     counting and the exact pair was repeated at commit `88c018086b21`.
   - Run cached and uncached VEP 115.2 `--everything` on all 50 1000 Genomes
     samples.
   - Use one warm-up and three measured repetitions, randomized within sample.
5. **Mechanism and complexity experiments** — pending
   - Measure controlled hit rates of 0%, 25%, 50%, 75%, 90%, 95%, and 100%.
   - Compare lightweight VEP, VEP `--everything`, a pinned plugin-rich clinical
     stack, and a lookup-only negative control.
   - Measure the real-input-size curve with matched WGS, Twist WES, and ACMG
     panel files. Use repeated short runs for the panel to estimate startup
     overhead and the break-even record count.
   - Run one paired annotation per HPRC R2 sample and report this graph-derived
     cohort separately as a generalization check.
6. **Analysis and manuscript** — pending
   - Treat samples as the statistical unit and report paired ratios, medians,
     IQRs, and bootstrap 95% confidence intervals.
   - Publish raw metrics, manifests, environment locks, plotting code, and figure
     source data.

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
