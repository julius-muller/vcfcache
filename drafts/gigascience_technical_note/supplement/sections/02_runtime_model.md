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
every model parameter from independently timed stages. Existing perfect-hit
diagnostics are described under fastVEP engineering pilots.
