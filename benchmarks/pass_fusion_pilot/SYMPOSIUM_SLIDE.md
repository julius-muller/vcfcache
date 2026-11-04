# Same project does not mean uniform performance

## Slide-ready wording

Even within one OpenStack project, VMs with the same advertised vCPU and RAM
allocation can expose different CPU models and deliver different performance.
Our Slurm workers exposed both Intel Broadwell and AMD EPYC-Rome CPUs; worker
assignment was associated with substantial differences in annotation runtimes.
For benchmarking and capacity planning, record the exposed CPU and storage
configuration and use same-node, counterbalanced comparisons—not resource
counts alone.

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
cache state, caused the difference. The new pinned-node pilot records these
properties explicitly; do not present optimistic modeled gains as measured ones.

## Official background

- [OpenStack Nova: CPU models](https://docs.openstack.org/nova/2026.1/admin/cpu-models.html)
  explains host-model, host-passthrough and configured guest CPU models.
- [OpenStack Nova: CPU topologies](https://docs.openstack.org/nova/latest/admin/cpu-topologies.html)
  describes shared/dedicated CPU allocation and topology policies. Availability
  depends on the provider; no infrastructure policy change is implied here.

Draft only: no slide, message or presentation has been published or sent.
