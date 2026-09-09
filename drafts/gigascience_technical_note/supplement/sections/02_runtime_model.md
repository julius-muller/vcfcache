## Runtime model and limiting cases

The main-text model separates non-annotation work in direct and cached paths
and starts the native annotator only if $M=N(1-f)>0$. Its general break-even
condition is

$$
A_1fN+(1-\mathbf{1}_{M>0})A_0>B_c(N)-B_d(N).
$$

Here $A_0$ is native startup, $A_1$ is average per-allele annotation time and
$B_c-B_d$ is additional processing cost. Perfect coverage removes startup
as well as per-allele work. With at least one miss and dominant annotation
cost, relative runtime approaches $1-f$ and speedup approaches $1/(1-f)$.
At $f=1$, the remaining $B_c(N)$ prevents infinite throughput.

The common $A_1$ is an approximation. If hits and misses differ in their
transcript overlap, variant types or supplementary-database matches, the
recoverable work is the sum of annotation costs over hits, not necessarily
$f$ times the total. Similarly, processing costs can depend on cache size,
output payload, storage and thread allocation as well as input count.
The curves in Figure 3 are explanatory approximations; they do not estimate
every model parameter from independently timed stages.

### Observed stage timings from an existing real-WGS run

The archived measured VEP run for SGDP genome `LP6005519-DNA_B04` with the
gnomAD AF ≥1% cache contains 4,876,494 input alleles and 342,007 misses
(92.99% hits). Its command durations provide the following breakdown without
rerunning the sample.

| Operation | Wall time (s) |
|---|---:|
| Multiallelic splitting and input filtering | 34.409 |
| Cache annotation lookup | 102.368 |
| Miss extraction | 44.556 |
| VEP on misses, including input/output conversion | 911.632 |
| Merge annotations | 52.822 |
| Final filtering and output | 59.591 |
| Other end-to-end overhead (difference) | 1.902 |
| End-to-end cached runtime | 1,207.280 |

Each command includes its own output processing, compression and inline index
creation where specified. The six sequential command durations sum to
1,205.378 s. The residual is calculated against the recorded end-to-end time,
not attributed to a separately measured component. Direct annotation took
11,467.473 s; the retained pair passed its VEP semantic gate and achieved
9.499-fold speedup. Thus real-sample stage timing already exists; only
separation of internal costs such as compression versus lookup would require
additional instrumentation.

Evidence: `figures/source_data/external_wgs_metrics.tsv`, VEP/sample/AF ≥1% row,
and the archived `workflow.log` under
`benchmarks/results/raw/campaigns/external-vep-light-full52-bbef41d7db8e/measured/tasks/task-23/runs/gnomad_af_0.01/pilot/LP6005519-DNA_B04/bbef41d7db8e/cached_r01/stats/LP6005519-DNA_B04_vcstats/`.
The final native campaign's 180 cached cells contain no zero-miss observation;
the perfect-hit case discussed in the model and early engineering pilot is a
separate controlled condition.
