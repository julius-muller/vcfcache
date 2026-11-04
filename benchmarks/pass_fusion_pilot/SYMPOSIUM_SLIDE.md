# Same project does not mean uniform performance

## Slide-ready wording

Even within one OpenStack project, VMs with the same advertised vCPU and RAM
allocation can expose different CPU models and deliver different performance.
In our controlled Slurm test, the same annotation workload allocated eight CPUs
took about 51 seconds on one worker versus 84–98 seconds on others—up to a
1.9-fold runtime difference.
For reproducible benchmarks and reliable capacity planning, record the exposed
CPU and storage configuration, and compare alternatives on the same worker
with balanced run order—not resource counts alone.

## Short takeaway

**Benchmark the instance, not just the flavor.**

## Speaker note for the cBioPortal infrastructure talk

This is a general infrastructure lesson from our annotation/Slurm benchmark,
not a measured cBioPortal performance comparison. Hardware heterogeneity can
affect capacity planning as well as experiments: a software version can appear
faster merely because its control ran on another worker. For cBioPortal, test
representative queries and imports on the intended deployment, and record the
instance configuration alongside the software configuration.

In a live check, `sl-w2` reported `Intel Core Processor (Broadwell, IBRS)` and
`sl-w5` reported `AMD EPYC-Rome Processor`. Both exposed 16 guest CPUs, local
`/mnt/data` on ext4 `/dev/sdb`, and a shared NFS `/results` mount. Guest labels
do not prove exact physical-hypervisor specifications. The historical timing
association alone does not establish that storage, rather than CPU, memory or
cache state, caused the difference. The controlled pilot used KPGP-00139,
vep-rs core, fresh direct annotation, eight-CPU Slurm affinity, and two reversed-
order technical repeats per worker. Node medians were 51.34 s on `sl-w5` and
83.71–98.03 s on the other five workers. These are repeated measurements of one
workload, not independent biological replicates or a causal CPU-model study.
Both direct and cached runs were faster on `sl-w5`; attributing the difference
exclusively to disk speed would be unjustified. Guest CPU models were AMD
EPYC-Rome on `sl-w5` and Intel Broadwell on the other workers.

The numbers are in `results-20261104-v2/conditions.tsv`; the frozen source is
`e72d6fe11fc2df302bea6dd6f7987d0b0fbb08d0`. This is an annotation/Slurm example
to motivate measuring cBioPortal imports and queries on their intended VMs,
not evidence of a measured 1.9-fold cBioPortal difference.

## Official background

- [OpenStack Nova: CPU models](https://docs.openstack.org/nova/2026.1/admin/cpu-models.html)
  explains host-model, host-passthrough and configured guest CPU models.
- [OpenStack Nova: CPU topologies](https://docs.openstack.org/nova/latest/admin/cpu-topologies.html)
  describes shared/dedicated CPU allocation and topology policies. Availability
  depends on the provider; no infrastructure policy change is implied here.

Draft only: no slide, message or presentation has been published or sent.
