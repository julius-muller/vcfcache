** WARNING: connection is not using a post-quantum key exchange algorithm.
** This session may be vulnerable to "store now, decrypt later" attacks.
** The server may need to be upgraded. See https://openssh.com/pq.html
# Slurm worker CPU mismatch — 9 November 2026

Read-only live inspection from `itcc-monitoring`; no VM, Slurm or cloud configuration was changed.

| Nodes | Guest-visible CPU model | vCPUs per node |
|---|---|---:|
| sl-w1, sl-w2, sl-w3, sl-w4, sl-w6 | Intel Core Processor (Broadwell, IBRS) | 16 |
| sl-w5 | AMD EPYC-Rome Processor | 16 |
| sl-head | AMD EPYC-Rome Processor | 8 |

All six workers use the same `de.NBI medium` flavor (32 GiB RAM), Ubuntu
26.04 image and kernel `7.0.0-38-generic`. The five Intel workers match on
guest vendor/model/family/stepping, instruction flags and CPU topology. Worker
5 differs in vendor, model, flags and reported cache configuration. These are
virtual CPU descriptions, not verified physical processor SKUs; matching
descriptions alone cannot establish equal physical host performance.

The seven nodes have distinct project-scoped OpenStack host identifiers. The
worker group uses soft anti-affinity: it prefers separate hosts but does not
constrain CPU vendor or generation. The flavor's only placement extra spec is
`aggregate_instance_extra_specs:compute=true`; there is no CPU restriction.
Only the `nova` availability zone is exposed. The project account receives
HTTP 403 when listing host aggregates or hypervisors, and actual compute-host
names are redacted from server details.

All six workers were created on 31 July 2026. Worker 5's visible action history
contains creation/volume attachment that day and a stop/start on 23 October;
it contains no reported migration. This does not establish which CPU model
it exposed previously or what CPU a newly scheduled VM would receive today.
There is no evidence that AMD has become the only/default placement target.

All workers were idle and the queue was empty during inspection. Slurm has
no CPU features configured, so it currently treats all six workers as eligible
for the same benchmark partition. Existing timings involving worker 5 should
be evaluated separately; this inspection does not establish a numerical bias.

## Recovery recommendation

1. For immediate runs, explicitly exclude `sl-w5` in the submission command
   (`sbatch --exclude=sl-w5 ...` or `srun --exclude=sl-w5 ...`). Existing
   submissions without that option remain eligible to use it.
2. Ask de.NBI for placement on a matching Intel host, preferably through a
   CPU-specific flavor/host aggregate or equivalent scheduling constraint.
   Confirm the physical CPU model and whether dedicated/pinned cores are
   available for reproducible benchmarks. Do not merely mask AMD as Intel in
   the guest: a virtual CPU label is not equivalent hardware.
3. Once placement is available, drain only worker 5 and recheck that it has
   no allocated jobs. Prefer a provider-assisted move or, where supported,
   shelve/offload and unshelve onto the approved host. Verify cloud permissions,
   target availability and root snapshot completion before relying on this
   route; it has not been exercised in this audit. Rescheduling can also
   select AMD again without a placement constraint.
4. If recreation is necessary, preserve its existing 500 GiB data volume and
   back up root configuration/software first. The volume has
   `delete_on_termination=false`, contains approximately 309 GiB of data, and
   must be reattached and mounted without formatting. `/results` is the
   controller's NFS export. Reconcile the authoritative Slurm Terraform state
   before a reviewed replacement plan: this checkout has the Slurm source
   but no local Slurm state or generated inventory. Do not apply it as a fresh
   stack or use the separate production-service Terraform state.
5. Recheck guest CPU details and provider placement; restore/verify software,
   mounts, Slurm registration and result access; run a representative
   performance check before returning the worker to benchmark service.
   Record hardware information with each benchmark run and recheck it after
   migrations/restarts. Keep head-node work outside measured worker execution.

Project quota is already **14/14 instances** (192/214 vCPUs allocated).
Consequently a parallel replacement/probe VM cannot currently be launched,
even though CPU/RAM quota has headroom. A temporary instance-quota increase
would permit a probe while retaining the original worker. Stopping or
shelving an instance does not remove it from the instance count.

Recreating all workers with the same unconstrained flavor could reproduce a
mixed cluster. It is unnecessary unless an alternative homogeneous platform
is deliberately selected and verified; retain all worker volumes and the
controller's data if that is later chosen. No replacement was attempted.

## Evidence and references

Protected live evidence:
`/home/ubuntu/slurm-audit-artifacts/20261109T110903Z/` (hardware, comparison
fields, worker identity/volume attachment, cloud inventory, quota and queue).

- [OpenStack CPU models](https://docs.openstack.org/nova/latest/admin/cpu-models.html)
- [Host aggregates and flavor constraints](https://docs.openstack.org/nova/latest/admin/aggregates)
- [Scheduling and optional same-host/different-host hints](https://docs.openstack.org/nova/latest/admin/scheduling.html)
- [Server lifecycle: shelve/offload/unshelve](https://docs.openstack.org/api-guide/compute/server_concepts.html)

OpenStack same-host/different-host scheduler hints are possible only when the
corresponding filters are enabled; their availability here is unverified.
Co-locating a replacement with an existing Intel worker can also introduce
resource contention and conflict with the intended worker separation, so it
is not the preferred benchmark fix.

## Outcome — sl-w5 replaced (same day)

sl-w5 was drained, its root configuration archived to
`sl-head:/mnt/data/backups/sl-w5-root-20261109/` (sha256 `ac06b230…`), its
data volume `b782955d-375d-456c-a7eb-cd886e2011b4` detached intact, and the
AMD VM deleted. In the freed slot, two throwaway `de.NBI medium` probes and
the new worker all landed on Intel Broadwell, on three different hosts. The
new sl-w5 (`b02cec31-7c8d-4c05-9a06-d21d4d651b71`, **10.0.0.229**, previously
10.0.0.159) was kept, the original volume was reattached unformatted (309 GiB
preserved), and the Slurm role was applied to the whole cluster. The only
change on the other nodes was the sl-w5 address in `/etc/hosts`/`slurm.conf`
plus daemon restarts; the `/mnt/data/vcfcache` mode reset by the role was
restored to the live 0775.

Post-replacement check (`openssl speed sha256`, 16 KiB blocks, via Slurm):

| Node | Kernel | CPU | 1 thread (kB/s) | 16 threads (kB/s) |
|---|---|---|---:|---:|
| sl-w1 | 7.0.0-38 | Broadwell | 324626 | 5159201 |
| sl-w2 | 7.0.0-38 | Broadwell | 324298 | 5198470 |
| sl-w3 | 7.0.0-38 | Broadwell | 313365 | 4845083 |
| sl-w4 | 7.0.0-38 | Broadwell | 322919 | 5033830 |
| sl-w5 | 7.0.0-38 | Broadwell | 321811 | 5161111 |
| sl-w6 | 7.0.0-38 | Broadwell | 325254 | 5164640 |

Residual differences and caveats:

- The cluster image `Ubuntu 26.04 LTS (2026-04-21)` was **deactivated by
  de.NBI on 2026-10-05**; sl-w5 was built from `Ubuntu 26.04 LTS (2026-09-27)`
  (`0f95cdb1-…`). It was pinned to the cluster's running kernel 7.0.0-38 via
  `/etc/default/grub.d/90-slurm-kernel.cfg`. Userland packages differ
  slightly in both directions (newer image base vs. security updates on the
  older workers). The other workers have 7.0.0-42 installed and will boot it
  on their next reboot; sl-w5 will not until the GRUB pin is removed.
- The Apptainer PPA no longer publishes the pinned `1.5.3-1~resolute`, and
  Launchpad returns 404 for that build. sl-w5's package was rebuilt from
  sl-head's installed copy with `dpkg-repack`; all 274 files match sl-w1's
  dpkg md5sums. It is on hold. The .deb remains in
  `sl-head:~/apptainer-repack/` (sha256 `56b31ad7…`). A fresh role run on a
  new node will fail at that task until the pin or the source is changed.
- The `ansible.builtin.command`/`shell` probes in the role now set
  `check_mode: false` so that `--check --diff` works.
- There is still no Terraform state for this stack. sl-w5 was created with
  the OpenStack CLI using the same name, flavor, keypair, security group,
  server group and metadata as `main.tf`. `cluster_image_name` in
  `variables.tf` now names an inactive image. Do not `terraform apply` this
  root without importing the live resources first.
- CPU placement is still unconstrained; a future rebuild/migration can land
  on AMD again. Recheck `lscpu` after any reprovisioning.

Evidence: `/home/ubuntu/slurm-audit-artifacts/20261109T110903Z/`
(`sl-w5-predelete.json`, `placement-probes.txt`, `post-replacement-cpucheck.txt`).

## Update — full package upgrade (same day)

All seven nodes received `apt-get update` + `full-upgrade` (phased
`software-properties` 0.120.1 installed explicitly so all nodes match) and were
rebooted, controller first. This supersedes the sl-w5 kernel pin and
repacked Apptainer described above: the GRUB pin was removed, the hold
released, and the role now pins **Apptainer 1.5.5-1~resolute**; its version
check follows the variable instead of a hard-coded "1.5.3".

State afterwards: every node runs kernel **7.0.0-43-generic**, Apptainer
1.5.5, Slurm 25.11.2-1ubuntu2 (unchanged, apt-pinned); no upgradable packages
remain; all data volumes and `/results` mounted; a full role run
is clean. The role's `/mnt/data/vcfcache` mode is now 0775 to match the live
nodes.

Repeated `openssl speed sha256` runs on the Intel workers vary by up to ~15%
single-threaded and ~8% at 16 threads between consecutive runs on the same node
(e.g. sl-w2 1T 274k–325k kB/s), consistent with shared-host contention rather
than a hardware class difference. Benchmarks should use repeated runs and
record per-run node/CPU/kernel. Evidence: `post-upgrade-cpucheck.txt`.
