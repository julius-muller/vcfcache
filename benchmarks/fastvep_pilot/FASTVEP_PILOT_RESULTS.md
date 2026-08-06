# Exploratory fastVEP pilot results

> Follow-up: lightweight statistics and a ten-database supplementary stress
> profile are now complete. The heavier pipeline reached 3.60x at 90% hits and
> 8.57x at 100% hits with complete output equality. See
> [HEAVY_FASTVEP_RESULTS.md](HEAVY_FASTVEP_RESULTS.md). The results below remain
> the frozen original core-pipeline pilot with the former full statistics cost.

## Conclusion

Current VCFcache is already useful with fastVEP, but it does not yet meet the
new fast-annotator product target. On one 4,329,621-variant WGS, VCFcache
reduced end-to-end runtime by 46% at 80% hits and 53% at 90% hits while
preserving the complete output. This corresponds to 1.86x and 2.11x speedups,
respectively. The target was 2x at 80% and 3x at 90%.

The 100%-hit diagnostic reached 2.84x. VCFcache therefore provides a material
benefit even around this fast native annotator, but fixed cache-processing and
statistics costs cap the achievable speedup. The native/streaming work should
remain a separate follow-up project until the cheapest improvements have been
evaluated.

fastVEP itself should currently be classified as **promising but incomplete**,
not as a drop-in replacement for VEP. It retained every smoke-test variant and
was highly concordant on most shared core fields, but transcript coverage
differed substantially, HGVSc concordance was 98.14%, and the advertised
`--symbol`, `--canonical`, and `--everything` switches had no effect in the
pinned version.

## Configuration

- Date: 2026-08-05.
- Host: `gvbrowse-preproc` in ITCCcloud_dev; 16 Rayon/bcftools threads.
- VCFcache: 0.5.2 from the existing benchmark VM checkout.
- fastVEP: 0.3.0 at commit
  `e47216cebe3abcd8dff722b7fb0ab1b19d4fcc80`.
- Transcript input: Ensembl release 115 GRCh38 GFF3 and matching FASTA.
- WGS: normalized, biallelic HG02079 GRCh38 input, 4,329,621 variants.
- Profile: fastVEP VCF output with HGVS. The CLI also received `--symbol`,
  `--canonical`, and `--everything`, which the flag smoke showed to be no-ops.
- No technical repeats. Cache construction and fastVEP transcript-cache
  construction were outside timed cells.

Both arms used the same input and `--skip-split-multiallelic`; the frozen input
contained neither multiallelic nor spanning-deletion records. Default VCFcache
statistics were included in every end-to-end time.

## WGS performance and output equality

| Condition | Observed hits | Misses | Wall time | Speedup | Workflow-only speedup | Exact output |
|---|---:|---:|---:|---:|---:|:---:|
| Direct | 0% | 4,329,621 | 389.2 s | 1.00x | 1.00x | reference |
| 80% target | 80.02% | 865,183 | 209.0 s | 1.86x | 2.16x | yes |
| 90% target | 90.02% | 432,201 | 184.6 s | 2.11x | 2.57x | yes |
| 100% target | 100.00% | 0 | 136.9 s | 2.84x | 4.64x | yes |

All three cached outputs matched the direct output across all 4,329,621 ordered
records after canonicalizing only INFO-tag and CSQ-entry order. CSQ headers were
equal, there were zero record mismatches, and every comparison produced the
same canonical SHA-256:
`15c5d795296ec1598237b907f1f3f1af0bf6bb70075ce00bfce66fb90fdacbf6`.

The end-to-end 2x/3x performance gate therefore **failed**, while the output
equality gate passed.

## Where time is spent

At 100% hits, no fastVEP annotation was executed. The 72.7-second VCFcache
workflow comprised approximately:

- 5.4 s input filtering;
- 24.5 s cache annotation lookup;
- 20.1 s miss selection;
- 0.9 s no-miss merge/copy;
- 21.6 s final filtering and BCF writing.

Default statistics and surrounding CLI work added another 64.2 seconds to the
100%-hit end-to-end time. Similar post-workflow overhead was 51.9–53.3 seconds
in the other cells. The internal workflow traces imply 2.16x, 2.57x, and 4.64x
speedups at 80%, 90%, and 100%, respectively, but these are diagnostic ratios,
not separately measured `--no-stats` runs.

The immediate optimization order is therefore:

1. measure a matched `--no-stats` or lightweight-statistics mode;
2. eliminate redundant full-file statistics scans for ordinary users;
3. reduce the lookup, miss-selection, and final-output passes;
4. only then decide whether a generic streaming engine, standalone fastVEP
   wrapper, or upstream/fork integration is warranted.

## fastVEP versus VEP smoke

The smoke used 8,769 variants from chr22:20–25 Mb. fastVEP and VEP retained all
8,769 variants. The comparison deliberately considered complete transcript
coverage before field concordance:

- fastVEP allele-transcript pairs: 158,964;
- VEP allele-transcript pairs: 72,095;
- shared pairs: 70,141;
- fastVEP-only pairs: 88,823;
- VEP-only pairs: 1,954.

Among shared pairs, Allele, gene/transcript identifiers, BIOTYPE, canonical
status, CDS/protein position, and several other core fields were identical.
Consequence matched for 99.996%, HGVSp for 99.987%, and HGVSc for 98.139%.
These numbers do not explain whether differences are fastVEP errors, VEP GFF
filtering behavior, or intentionally different transcript selection. That
requires a focused validation before any replacement claim.

Adding `--symbol --canonical --everything` to the fastVEP `--hgvs` command
produced identical CSQ headers and a byte-identical VCF body. This agrees with
the static audit showing those values are currently parsed but not forwarded
by the CLI.

## Interpretation boundary

This was an exploratory engineering pilot on one machine and one WGS. It is
sufficient to choose the next engineering work, but not to make a publication
claim about fastVEP correctness or a hardware-independent performance claim.
The larger ANNOVAR/SnpEff/fastVEP campaign remains gated on this decision.

Machine-readable source tables are under `source_data/`. Raw logs, BCFs,
resource reports, equality reports, and the environment manifest remain at
`/mnt/data/vcfcache_benchmarks/fastvep_pilot` on ITCCcloud_dev.
