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

![](supplement/figures/supplementary_figure3_matched_assays_v3.png){width=100%}

**Supplementary Figure 3: Sample-level cache coverage and performance across
matched assay scales.** Twelve independent GRCh38 genomes were represented as
ACMG SF v3.3 Panel, Twist Human Core Exome WES and WGS inputs. Colour and shape
identify the provider cohort; black horizontal marks denote medians. (A) The
fraction of variants found in the gnomAD AF ≥1% cache. Cache membership is
annotator-independent and is therefore shown once per sample and assay. (B)
End-to-end speedup for VEP 115.2 and fastVEP 0.3.0. The dashed line denotes
equal direct and cached runtime; values below it identify cells in which fixed
cache-workflow costs exceeded annotation work saved. Each sample-tool-assay
combination was measured once, and every cached output passed its
annotator-specific semantic comparator. ACMG, American College of Medical
Genetics and Genomics; AF, allele frequency; KPGP, Korean Personal Genome
Project; SGDP, Simons Genome Diversity Project; VEP, Variant Effect Predictor;
WES, whole-exome sequencing; WGS, whole-genome sequencing.
