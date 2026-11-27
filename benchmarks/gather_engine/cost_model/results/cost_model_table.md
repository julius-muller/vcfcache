| Annotator | Sources | Cache | Miss-step blocks | Direct (s) | Cached (s) | Misses | Speed-up (range) | Ideal | Share of ideal | Cost per miss |
|---|---|---|---|---|---|---|---|---|---|---|
| fastVEP 0.5.0 | 11 | gnomAD AF>=1% | 65 kb | 681 | 336 | 8.2% | 2.03x (2.02-2.06) | 5.8x | 35% | 4.86x |
| fastVEP 0.5.0 | 11 | gnomAD AF>=1% | 4 kb | 681 | 270 | 8.2% | 2.52x (2.50-2.54) | 5.8x | 43% | 3.70x |
| fastVEP 0.5.0 | 11 | gnomAD AF>=1% | 1 kb | 681 | 209 | 8.2% | 3.26x (3.26-3.29) | 5.8x | 56% | 2.61x |
| fastVEP 0.5.0 | 11 | cohort | 65 kb | 681 | 397 | 2.9% | 1.73x (1.72-1.74) | 8.3x | 21% | 16.63x |
| fastVEP 0.5.0 | 11 | cohort | 4 kb | 681 | 201 | 2.9% | 3.40x (3.24-3.41) | 8.3x | 41% | 6.98x |
| fastVEP 0.5.0 | 11 | cohort | 1 kb | 681 | 161 | 2.9% | 4.23x (4.06-4.26) | 8.3x | 51% | 4.99x |
| fastVEP 0.5.0 | 11 | rerun (all hits) | - | 681 | 62 | 0.0% | 11.16x (10.81-11.21) | 11.2x | 100% | - |
| fastVEP 0.5.0 | 15 | gnomAD AF>=1% | 65 kb | 998 | 537 | 8.2% | 1.80x (1.73-1.86) | 6.5x | 28% | 5.80x |
| fastVEP 0.5.0 | 15 | gnomAD AF>=1% | 4 kb | 998 | 414 | 8.2% | 2.38x (2.27-2.41) | 6.5x | 36% | 4.29x |
| fastVEP 0.5.0 | 15 | gnomAD AF>=1% | 1 kb | 998 | 288 | 8.2% | 3.37x (3.32-3.46) | 6.5x | 52% | 2.74x |
| fastVEP 0.5.0 | 15 | cohort | 65 kb | 998 | 587 | 2.9% | 1.64x (1.56-1.71) | 10.0x | 16% | 17.66x |
| fastVEP 0.5.0 | 15 | cohort | 4 kb | 998 | 279 | 2.9% | 3.48x (3.35-3.57) | 10.0x | 35% | 7.39x |
| fastVEP 0.5.0 | 15 | cohort | 1 kb | 998 | 229 | 2.9% | 4.30x (4.29-4.36) | 10.0x | 43% | 5.39x |
| fastVEP 0.5.0 | 15 | rerun (all hits) | - | 998 | 67 | 0.0% | 13.98x (13.76-14.98) | 14.0x | 100% | - |
| vep-rs 0.5.0 | 11 | gnomAD AF>=1% | - | 2937 | 454 | 8.2% | 6.40x (6.26-6.61) | 10.1x | 64% | 1.65x |
| vep-rs 0.5.0 | 11 | cohort | - | 2937 | 250 | 2.9% | 11.80x (11.46-12.26) | 20.8x | 57% | 2.23x |
| vep-rs 0.5.0 | 11 | rerun (all hits) | - | 2937 | 51 | 0.0% | 56.22x (52.22-58.39) | 56.2x | 100% | - |
| vep-rs 0.5.0 | 15 | gnomAD AF>=1% | - | 6391 | 957 | 8.2% | 6.56x (6.47-6.68) | 11.0x | 60% | 1.72x |
| vep-rs 0.5.0 | 15 | cohort | - | 6391 | 463 | 2.9% | 13.09x (12.83-13.98) | 26.0x | 50% | 2.20x |
| vep-rs 0.5.0 | 15 | rerun (all hits) | - | 6391 | 56 | 0.0% | 108.76x (106.48-114.76) | 108.8x | 100% | - |
| Ensembl VEP 116.3 | 11 | gnomAD AF>=1% | - | 13211 | 1261 | 8.2% | 10.54x (10.41-10.71) | 11.6x | 91% | 1.10x |
| Ensembl VEP 116.3 | 11 | cohort | - | 13211 | 536 | 2.9% | 24.52x (23.58-24.71) | 30.1x | 81% | 1.25x |
| Ensembl VEP 116.3 | 11 | rerun (all hits) | - | 13211 | 56 | 0.0% | 238.45x (235.59-240.87) | 238.5x | 100% | - |
| SnpEff 5.4c | 8 | gnomAD AF>=1% | - | 5818 | 830 | 8.2% | 6.92x (5.86-7.23) | 10.0x | 69% | 1.50x |
| SnpEff 5.4c | 8 | cohort | - | 5818 | 574 | 2.9% | 10.02x (9.88-10.83) | 20.6x | 49% | 2.66x |
| SnpEff 5.4c | 8 | rerun (all hits) | - | 5818 | 108 | 0.0% | 54.95x (51.29-56.73) | 55.0x | 100% | - |

Whole genomes (KPGP-00006, -00086, -00128; 4.6-4.8 M variants), 3 repeats each, median over repeats then over samples; outputs identical to direct annotation in every run, except that Ensembl VEP fills HGNC_ID differently in 164-277 records per genome (at most 0.006%; known VEP bug Ensembl/ensembl-vep#1959, no other field differs) and that SnpEff's cached output lists the same INFO tags and values in another order in 6-7% of records. SnpEff runs the 8 VCF sources of the Ensembl VEP stack through SnpSift. Direct = the annotator alone as advised (fastVEP sources at 65 kb blocks). Ideal = D / (R + m D): the speed-up if every miss cost what one variant costs in direct annotation, with R the rerun time (vcfcache's own cost). Cost per miss = (C - R) / (m D), relative to that.
