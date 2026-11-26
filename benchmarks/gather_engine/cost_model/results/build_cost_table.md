| Annotator | Cache build (h) | Genome: direct / cached (s) | Break-even (genomes) | Exome: direct / cached (s) | Break-even (exomes) |
|---|---|---|---|---|---|
| fastVEP 0.5.0 | 0.21 | 681 / 209 | 1.6 | 242 / 64 | 4.3 |
| vep-rs 0.5.0 | 1.54 | 2937 / 454 | 2.2 | 162 / 53 | 50.7 |
| Ensembl VEP 116.3 | 10.12 | 13211 / 1261 | 3.0 | 510 / 189 | 113.5 |

Cache build = annotating the 18.9 M-record gnomAD AF>=1% blueprint with the 11-source stack on one node (16 threads; Ensembl VEP: the sum of 12 parts, one node each). Break-even = build / (direct - cached): the number of samples after which building the cache for one's own stack has paid for itself. With a published cache there is no build. fastVEP's 1 kb miss-step sources are a one-time build of their own, not included.
