## Native clinical annotation campaign

The final campaign extended the assay-scale design to SnpEff and to clinical
fastVEP recipes on both assemblies. Recipes were compared only with their own
direct baseline because VEP `--everything`, fastVEP and SnpEff expose different
annotation content. The clinical SnpEff recipe comprised SnpEff ANN consequence
prediction followed by SnpSift ClinVar and gnomAD exome AF/AC/AN annotation.
The core-to-clinical SnpEff comparison on GRCh37 held samples and cache
membership fixed and therefore isolates additional per-variant annotation work.

**Supplementary Table 1: End-to-end speedup in the final native annotation
campaign.** Values are medians across concordant sample pairs. The four
warning-only SnpEff exclusions described above reduce two GRCh37 WES and two
GRCh37 WGS rows from 12 to 11 genomes. AF, allele frequency; WES, whole-exome
sequencing; WGS, whole-genome sequencing.

| Tool | Recipe | Assembly | Assay | Genomes | Median hit rate | Median speedup |
|---|---|---|---|---:|---:|---:|
| fastVEP | clinical | GRCh38 | Panel | 12 | 90.4% | 0.46-fold |
| fastVEP | clinical | GRCh38 | WES | 12 | 93.2% | 0.93-fold |
| fastVEP | clinical | GRCh38 | WGS | 12 | 92.6% | 2.54-fold |
| fastVEP | clinical | GRCh37 | Panel | 12 | 90.7% | 0.41-fold |
| fastVEP | clinical | GRCh37 | WES | 12 | 94.6% | 1.02-fold |
| fastVEP | clinical | GRCh37 | WGS | 12 | 95.8% | 1.70-fold |
| SnpEff | core | GRCh37 | Panel | 12 | 90.7% | 1.07-fold |
| SnpEff | core | GRCh37 | WES | 11 | 94.7% | 1.16-fold |
| SnpEff | core | GRCh37 | WGS | 11 | 95.9% | 4.46-fold |
| SnpEff | clinical | GRCh38 | Panel | 12 | 90.4% | 1.06-fold |
| SnpEff | clinical | GRCh38 | WES | 12 | 93.2% | 2.66-fold |
| SnpEff | clinical | GRCh38 | WGS | 12 | 92.6% | 10.78-fold |
| SnpEff | clinical | GRCh37 | Panel | 12 | 90.7% | 1.03-fold |
| SnpEff | clinical | GRCh37 | WES | 11 | 94.7% | 3.62-fold |
| SnpEff | clinical | GRCh37 | WGS | 11 | 95.9% | 12.31-fold |

Panels were at or below practical break-even. Every whole-genome recipe
benefited, while clinical SnpEff also benefited materially for exomes. On the
same GRCh37 inputs, adding clinical lookups increased median speedup from 1.16-
to 3.62-fold for WES and from 4.46- to 12.31-fold for WGS. This supports the
runtime model's central prediction: cache reuse becomes more valuable as the
recoverable work performed for each variant increases.
