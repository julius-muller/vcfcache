# Route-mechanism follow-up, job 1255

Copied read-only from the Slurm controller during adaptive-routing development.
The TSVs and original runner are raw experiment evidence, not maintained source.
Use worker measurements, not head-node preflights, for within-experiment contrasts.

This follow-up investigates cache volume, contig ordering, thread allocation and
group imbalance. Its input cells must not be treated as biological coverage or
publication-cohort results. The adaptive-routing validation is a separate
experiment under `../vcfcache-adaptive-routing-20261023/`.
