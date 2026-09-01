# Supplementary Figures

![](supplement/figures/supplementary_figure1_external_wgs_v2.png){width=100%}

**Supplementary Figure 1: Cache coverage and wall time returned across 52
independent WGS samples.** Each symbol represents one evaluation genome; black
horizontal marks denote medians. Colour and shape both identify the provider
cohort. (A) Proportion of input variants found in each cache strategy. Cache
membership is annotator-independent, so each genome is shown once. (B)
Percentage of end-to-end wall time saved with VEP 115.2 and fastVEP 0.3.0.
Public caches were downloaded from Zenodo. Cohort caches were constructed from
three provider-specific genomes disjoint from evaluation. AF, allele
frequency; KPGP, Korean Personal Genome Project; PGP, Harvard Personal Genome
Project; SGDP, Simons Genome Diversity Project; VEP, Variant Effect Predictor;
WGS, whole-genome sequencing.

![](supplement/figures/supplementary_figure2_fastvep_cpu_complexity.png){width=100%}

**Supplementary Figure 2: CPU availability does not eliminate annotation-reuse
benefit in fastVEP.** The same held-out 4.80-million-variant KPGP genome was
processed using the core consequence recipe and an engineering stress recipe
containing ten supplementary databases. A controlled cache constructed only
from the corresponding direct annotations fixed cache coverage at approximately
90%; it was used to isolate the interaction between annotation work and CPU
availability, not to estimate biological cache coverage. Process-wide CPU
affinity was restricted to 1, 10 or 32 CPUs for the core recipe and to 10 or
32 CPUs for the dense recipe, enabling direct comparison of dense fastVEP at
10 CPUs with core fastVEP at one CPU. (A) End-to-end runtime for direct
fastVEP and VCFcache. (B) Speedup over direct fastVEP at the same CPU limit.
Each retained cell was run once and each cached output passed strict
complete-record comparison with its direct counterpart. KPGP, Korean Personal
Genome Project; VCF, variant call format; VCFcache, cached annotation workflow;
WGS, whole-genome sequencing.

![](supplement/figures/v4/supplementary_figure3_annotation_load_fidelity_v4.png){width=100%}

**Supplementary Figure 3: Annotation load increases reuse benefit without a
cache-attributable change in output.** (A) Paired SnpEff speedups for core and
clinical recipes on the same GRCh37 PGP samples. Lines connect identical
samples; the dashed line marks equal runtime. The documented warning-only
sample was excluded before inference for WES and WGS. (B) Of 180 direct/cached
cells in the final clinical campaign, 176 were eligible and semantically
concordant. Four SnpEff cells were excluded because direct SnpEff itself changed
output with input order; these are not counted as evidence for or against
VCFcache. No eligible pair contained a cache-attributable mismatch. PGP,
Personal Genome Project; WES, whole-exome sequencing; WGS, whole-genome
sequencing.
