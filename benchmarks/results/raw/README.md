# Raw evidence

Results and provenance only: metrics TSVs, run manifests, campaign definitions,
sample and interval manifests, cache-build timings, and the cluster/tool
configuration the numbers were produced against.

**Tool logs and per-task retry debris are deliberately not here.** They live on
the Slurm cluster under `/mnt/data/slurm-results` (`benchmark_archive/`,
`*/logs/`, `campaigns/*/attempts/`) — roughly 22,000 files that duplicate in git
what the cluster already stores. What each run did, and which job produced it, is
recorded in `../RUN_LOG.md` instead.

Nothing was deleted: the cluster copy is the archive of record.
