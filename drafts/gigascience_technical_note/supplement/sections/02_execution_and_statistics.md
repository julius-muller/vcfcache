## Benchmark execution and statistical analysis

The Slurm cluster comprised one head node and six workers with identical
software and worker-local benchmark inputs, caches and annotator assets. A
shared results export was mounted read-only or append-only as appropriate.
Each sample task executed its direct and cached conditions in separate Slurm
job steps. This prevented an earlier condition's cgroup memory peak from being
reported for a later condition.

The final campaigns used one observation for each unique
sample/annotator/cache combination. Samples provide biological and pipeline
heterogeneity; repeated executions of the same technical cell were not treated
as independent observations. A small repeatability control showed sub-percent
differences and was not pooled with the 52-genome analysis.

The external semantic comparator ran after timing and was excluded from all
reported wall times. No observation was calibrated or rescaled.

For the short matched Panel/WES extension, one separate KPGP-00319 assay subset
was processed before either timed condition in every task. This untimed step
loaded annotator, transcript and cache-lookup assets consistently; it was not
an evaluation-sample replicate and was excluded from all summaries. The need
for this control was established before retaining extension timings, when a
discarded pilot showed that cold-state condition order dominated short cells.

Medians summarize sample-level speedups. Stratified bootstrap intervals used
10,000 resamples within KPGP, SGDP and PGP while preserving a sampled genome's
paired conditions. This respects the deliberately heterogeneous 20/20/12
cohort allocation.

Relative speedup was interpreted within each assembly and annotator recipe.
In particular, the GRCh37 PGP configuration had a substantially faster direct
fastVEP baseline than the two GRCh38 configurations, whereas parsing, lookup
and output construction did not decrease in the same proportion. Its smaller
relative fastVEP speedup therefore reflects the larger fixed-cost fraction and
must not be interpreted as evidence of poorer biological cache coverage.
