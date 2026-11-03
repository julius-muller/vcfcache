# Figure captions

**Figure 1. Cold application-start caching performance in twelve WES inputs per
recipe.** (A) Each point is an independent input, coloured and shaped by cohort
(KPGP, Korean Personal Genome Project; SGDP, Simons Genome Diversity Project).
Values greater than one indicate faster cached annotation relative to that
input's own direct CLI baseline. The dashed reference marks break-even. Four
globally defined cached configurations are shown; no per-sample best choice is
made. (B) Mean within-input fractions of cached elapsed time attributable to
cache lookup, the miss-recipe invocation, and other CLI/I/O work. Recipe time
includes loading and recipe-owned processing, not just consequence calling.
Each invocation uses the same eight-CPU allocation and a separate Slurm step.
The first three inputs per recipe were screened direct-first; remaining inputs
used counterbalanced primary-arm order. Additional configurations reuse the
measured direct baselines, not newly retimed baselines. Processes start fresh;
filesystem caches are not flushed. These engineering measurements are not
technical-repeat averages or publication confidence intervals. WES, whole-exome
sequencing; CLI, command-line interface; I/O, input/output.

**Figure 2. Paired effect of lookup-worker allocation at a fixed eight-CPU
budget.** Lines connect the same independent WES input under four lookup workers
plus four shared cache-codec threads and eight lookup workers without an extra
cache-codec pool. Both configurations retain the input, not the whole cache,
and use the same recipe, cache, and output checks. Lower elapsed time is better.
The comparison changes both lookup concurrency and decoder mode; it does not
establish linear worker scaling or isolate either mechanism. All observations,
including regressions, are shown. Colours and symbols indicate cohorts as in
Figure 1. No assay-specific or annotator-specific dispatch policy was introduced.

Vector PDF/SVG and 600-dpi PNG exports are supplied. The two components of
Figure 1 are also exported separately in `figures-v3/panels/`. No whole-cache
warm-throughput figure is presented because the sessions did not pass the
capacity/validation gate.
