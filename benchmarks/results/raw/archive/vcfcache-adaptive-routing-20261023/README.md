# Adaptive-routing engineering validation, 2026-10-23

New experiment; does not replace publication campaign data or existing figures.

- Controls: jobs 1257 (8 CPUs, three diagnostic repetitions), 1258 (16 CPUs,
  one repetition), 1259 (1 CPU, one repetition).
- Real held-out SGDP WGS: jobs 1261 (`LP6005519-DNA_A06`, median-hit example)
  and 1262 (`LP6005592-DNA_C05`, lower-hit example).
- Controls use the existing VEP 115.2 everything gnomAD AF>=1% BCF cache.
  Real WGS uses the existing locally built clinical fastVEP cache derived from
  the bundled GRCh38 AF>=1% blueprint; it is not a bundled fastVEP download.
- Every timed condition uses light statistics. Complete-record/header audits
  are separate. There are no annotation-field exceptions.
- Controls are drawn from the cache and are **not cohort hit-rate findings**.
- Pilot A includes Slurm launch time in `measurements.json:wall_seconds`.
  Use each condition's `resources.json:wall_seconds` for CLI-only comparisons.
  Pilot B records CLI-only time directly, with orchestration time separately.
- Conditions have separate Slurm steps, affinity binding and cgroup snapshots.
  `largest_child_rss_kib` is not an aggregate memory peak; use `memory.peak`.
- Native CSI footprint estimates are not measured cache I/O.
- Source snapshots and executable hashes identify the measured deployments;
  they were frozen before the final documentation/telemetry-only edits.

Large public VCF/BCF artifacts and derived stores stay on the workers/controller,
not in Git. Small logs, statistics, comparisons and source snapshots are retained
here. Existing data are never rewritten by the summarizer.
