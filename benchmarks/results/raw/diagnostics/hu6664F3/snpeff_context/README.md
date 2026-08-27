# Minimal SnpEff 5.4c processing-order reproducer

This directory contains direct-SnpEff-only evidence for the sole warning-level
difference found in the final SnpEff VCFcache campaign. VCFcache is not invoked
by the reproducer.

## Result

Using SnpEff 5.4c and database `GRCh37.75`:

| input | target warning |
|---|---|
| `chr2:11273459 C>T` alone | absent |
| `chr2:11273458 A>AC`, then `chr2:11273459 C>T` | present |
| the same two records, target first | absent |

The warning is `WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS` on two `C2orf50`
transcript consequences. All other target `ANN` values are unchanged. The
sorted two-record input and its output are
`minimal_inputs/adjacent_insertion_then_target.bcf` and
`minimal_outputs/adjacent_insertion_then_target.bcf`. The deliberately reversed
input and its output are `minimal_inputs/target_then_adjacent_insertion.vcf` and
`minimal_outputs/target_then_adjacent_insertion.bcf`.

`minimal_matrix.tsv` contains the exact target `ANN` strings. The broader
`context_matrix.tsv` shows the same switch between a one-record blueprint input
and its ±1 kb neighborhood. Sample-derived contexts remain warning-free from
one record through the complete chromosome.

The three plain-VCF files under `issue_inputs/` are the compact public-ticket
reproducer. Job 1241 reran those files without VCFcache; `verification.tsv`
confirms that the insertion-first input alone carries the warning and that the
two warning-free target annotations have identical SHA-256 digests.

## Direct command shape

Each input was passed directly through the standard SnpEff invocation:

```bash
bcftools view INPUT -Ov \
  | java -Xmx8g -jar snpEff.jar -noLog -noStats -nodownload \
      -c snpEff.config -dataDir DATA GRCh37.75 - \
  | bcftools view -Ob -o OUTPUT
```

The tested JAR SHA-256 is
`5e8f75cbf908a33c6fb2e65c81e66fe31236cb21bb0195541c4703d8202c22b3`.
The complete installed toolchain and database hashes are recorded under
`benchmarks/reproduction/final_clinical_campaign/toolchain/`.

## Interpretation

SnpEff does not intentionally use the preceding variant's INFO fields for this
warning. Its source attaches the result of a transcript sanity check, but the
minimal inputs demonstrate that the emitted diagnostic metadata depends on
processing order. This is suitable for an upstream SnpEff issue. It is not
evidence of a VCFcache merge or lookup defect.

All files are frozen by `SHA256SUMS`. Job 1236 produced the context matrix; job
1238 produced the minimal reproducer; and job 1241 independently verified the
compact issue inputs. Job 1237 failed before annotation because of a
shell-script initialization error and produced no scientific result.
