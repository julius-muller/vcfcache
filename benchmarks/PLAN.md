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
3. **Benchmark runner** — implemented; pilot pending
   - Pin the VCFcache commit, container, bcftools, VEP 115.2 cache, CPUs, RAM,
     scratch storage, and thread counts.
   - Capture cgroup CPU, memory peak, I/O, wall time, variant counts, hit rate,
     cache size, and output concordance.
4. **Headline experiment** — pending pilot
   - Run cached and uncached VEP 115.2 `--everything` on all 50 1000 Genomes
     samples.
   - Use one warm-up and three measured repetitions, randomized within sample.
5. **Mechanism and complexity experiments** — pending
   - Measure controlled hit rates of 0%, 25%, 50%, 75%, 90%, 95%, and 100%.
   - Compare lightweight VEP, VEP `--everything`, a pinned plugin-rich clinical
     stack, and a lookup-only negative control.
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

Data preparation runs on the current 32-vCPU, 62-GiB VM. The repeated annotation
matrix should move to Slurm after one cached/uncached pilot determines per-job
CPU, RAM, scratch, and wall-time requirements.
