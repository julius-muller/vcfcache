# Dense fastVEP supplementary-annotation pilot

## Outcome

VCFcache remained strongly useful when fastVEP was configured with ten
supplementary databases. On the same 4,329,621-variant HG02079 WGS, direct
annotation took 707.4 seconds, while 90% and 100% cache hits took 196.2 and
82.6 seconds: **3.60x** and **8.57x** speedups. Both cached outputs matched the
direct output across all ordered records with zero mismatches and canonical
SHA-256
`ee5ac8103bff16110110aee14a54df4468db1a30ca93d15cfb3e94623899a7e4`.

| Condition | Misses | Wall time | Speedup | Complete output equality |
|---|---:|---:|---:|:---:|
| Direct | 4,329,621 | 707.4 s | 1.00x | reference |
| 90% hits | 432,201 | 196.2 s | 3.60x | yes |
| 100% hits | 0 | 82.6 s | 8.57x | yes |

These cells used the new `--statistics light` default. Cache construction and
the streaming equality check were outside timed cells. There were no technical
repeats.

## Stress configuration and boundary

The pilot loaded six dense `custom_vcf` fastSA databases, three dense
population-style databases built through fastVEP's gnomAD, 1000 Genomes, and
TOPMed parsers, and a real GRCh38 ClinVar database. The personal WGS was reused
as the source of dense deterministic records. This is deliberately an
upper-bound engineering stress configuration, not a biological population
reference or a manuscript estimate of ordinary supplementary-database density.

fastVEP 0.3.0 computed all six arbitrary custom sources and exposed them in
JSON, but did not project their arbitrary names into VCF INFO. The three
population sources and ClinVar did produce `FV_GNOMAD`, `FV_1KG`, `FV_TOPMED`,
and `FV_CLINVAR`, so the run stressed both supplementary lookup and VCF output
construction. This missing custom-to-VCF projection is a fastVEP output
limitation to track separately.

Direct heavy annotation took about 2.1 times the earlier 337.3-second core
fastVEP workflow. The additional burden is therefore real fastSA lookup and
output work, not an artificial delay or ineffective CLI flag.

## Clean 100%-hit micro-profile

The production-like commands measured 4.0 seconds for input filtering and BCF
writing, 22.2 seconds for cache lookup plus compressed output, 20.0 seconds for
the guaranteed-empty miss scan, and 18.3 seconds for final filtering plus
compressed output. Process startup was negligible: medians across 25 launches
were 1.7 ms for `/bin/true`, 11.9 ms for bcftools, and 3.0 ms for fastVEP.

Building a fresh index after writing either 523 MB output took about 15 seconds,
but the normal workflow uses bcftools `-W` to build the index inline, overlapping
that work with output writing. Null-output controls showed that cache matching
adds approximately 25.9 seconds over a plain BCF decode. Compression is
multithreaded and overlapped with record processing, so the measured wall-time
deltas are diagnostic rather than additive.

The remaining optimization target is consequently clear: remove the separate
miss-selection and final-output full-file passes, ideally by detecting the
all-hit case during lookup and assembling the final output in the same stream.
Eliminating subprocess startup would have no measurable effect.

Raw logs, BCFs, resource reports, source manifest, smoke gate, and equality
reports are stored at
`/mnt/data/vcfcache_benchmarks/fastvep_pilot/{heavy,profile_100_hit}` on
ITCCcloud_dev. Machine-readable summaries are under `source_data/` here.
