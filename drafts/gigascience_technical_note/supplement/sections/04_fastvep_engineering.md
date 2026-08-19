# Supplementary Results

## fastVEP engineering pilots

An exploratory pilot on one 4,329,621-variant WGS profiled the fixed work that
remains when cache hits approach 100%. With the original core fastVEP recipe,
VCFcache achieved 1.86-fold, 2.11-fold and 2.84-fold speedups at approximately
80%, 90% and 100% hits. Complete output comparison passed after canonicalizing
only INFO-tag and CSQ-entry order. This engineering pilot is not combined with
the primary external-WGS campaign.

A dense supplementary-database stress configuration used six custom fastSA
sources, three population-style sources and ClinVar. It was rerun on the
held-out KPGP-00319 genome under process-wide CPU affinity. At 10 CPUs, direct
and 90%-cached annotation took 706.2 and 230.0 seconds (3.07-fold speedup); at
32 CPUs, they took 691.4 and 205.4 seconds (3.37-fold). Increasing the limit
from 10 to 32 CPUs therefore changed direct runtime by only 2.1%, while
caching returned approximately eight minutes per genome at both limits. All
four retained outputs passed complete-record equality. This configuration
demonstrates that cache reuse remains valuable after fastVEP saturates its CPU
scaling, but its database density is deliberately an upper-bound stress test.

The matched core recipe took 868.4, 463.8 and 379.1 seconds directly at 1, 10
and 32 CPUs, and 576.5, 180.3 and 148.8 seconds with the controlled cache. The
corresponding speedups were 1.51-fold, 2.57-fold and 2.55-fold. Thus ten CPUs
made the dense direct recipe faster than the one-CPU core recipe (706.2 versus
868.4 seconds), but did not make complexity free: at the same 10-CPU allocation
the dense recipe was 52.2% slower than the core recipe. VCFcache returned 7.9
minutes from the dense 10-CPU run and increased its speedup from 2.57-fold for
the core recipe to 3.07-fold for the dense recipe. Additional CPUs and reuse
therefore addressed different components of runtime.

Profiling at 100% hits attributed material wall time to cache lookup, the
guaranteed-empty miss scan and final output construction; subprocess startup
was negligible. In the affinity-controlled experiment, both direct fastVEP and
VCFcache were restricted using `taskset`; the plotted axis therefore denotes
CPUs available to the process rather than only a configured thread count. An
exploratory dense one-CPU cell was stopped after sustained, impractically low
throughput and was not treated as a completed observation. It was unnecessary
for the prespecified comparison between the dense recipe at 10 CPUs and the
core recipe at one CPU.

The pilots motivate a future streaming engine that performs lookup, miss
selection and output construction with fewer complete-file passes. They do not
alter the current publication claim that VCFcache already reduces end-to-end
runtime around fastVEP.
