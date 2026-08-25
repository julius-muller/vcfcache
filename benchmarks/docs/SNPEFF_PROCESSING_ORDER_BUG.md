# SnpEff 5.4c processing-order-dependent transcript warning

## Summary

During strict direct-versus-cached validation, SnpEff 5.4c produced one
different diagnostic token on `chr2:11273459 C>T`. Direct SnpEff-only tests
reduced the behavior to two records and proved that it depends on processing
order. VCFcache is not involved in the reproducer.

With database `GRCh37.75`:

| direct SnpEff input | warning on `chr2:11273459 C>T` |
|---|---|
| target alone | absent |
| `chr2:11273458 A>AC`, then target | present |
| target, then `chr2:11273458 A>AC` | absent |

The order-dependent token is
`WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS`. It appears on the `C2orf50`
consequences for transcripts `ENST00000381585` and `ENST00000405022`. All other
target `ANN` content, including consequence, impact, gene, transcript and HGVS,
is unchanged.

The exact inputs, outputs, matrices and checksums are under
`benchmarks/results/raw/diagnostics/hu6664F3/snpeff_context/`. The compact
ticket inputs are in its `issue_inputs/` directory. SnpEff was invoked directly;
no cache lookup, VCFcache merge or VCFcache process was present.

## Environment

- SnpEff: 5.4c, build date 2026-02-23
- Database: `GRCh37.75`
- Java: Eclipse Temurin OpenJDK 21.0.12.1
- OS: Linux x86_64, kernel 7.0.0-14-generic
- SnpEff JAR SHA-256:
  `5e8f75cbf908a33c6fb2e65c81e66fe31236cb21bb0195541c4703d8202c22b3`

The direct command shape was:

```bash
java -Xmx8g -jar snpEff.jar \
  -noLog -noStats -nodownload \
  -c snpEff.config -dataDir /path/to/data \
  GRCh37.75 input.vcf > output.vcf
```

## Interpretation

SnpEff source attaches the result of `Transcript.sanityCheck()` to a variant
effect. The check describes the reference transcript protein and does not use
another input record's INFO annotations as evidence. Nevertheless, the minimal
pair proves that processing the adjacent insertion first changes the warning
later emitted for the target. This is consistent with unintended mutable or
cached transcript state inside SnpEff rather than intentional multi-variant
annotation.

The original VCFcache benchmark conservatively excludes the four cells in which
this raw warning difference occurred. Consequences, HGVS, ClinVar, gnomAD,
genotypes and record counts agree. VCFcache has not been changed to suppress or
rewrite the warning.

## Copy-paste GitHub issue

### Title

```text
SnpEff 5.4c: WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS depends on preceding VCF record order
```

### Body

````markdown
## Description

SnpEff 5.4c emits `WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS` for
`chr2:11273459 C>T` only when a nearby insertion is processed first. The target
alone is warning-free, and reversing the same two records removes the warning.
This was reproduced with direct SnpEff; no workflow cache or post-processing
tool is involved.

## Minimal inputs

All inputs use this header:

```vcf
##fileformat=VCFv4.2
##contig=<ID=chr2,length=243199373>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
```

Target alone:

```vcf
chr2	11273459	.	C	T	.	PASS	.
```

Insertion followed by target:

```vcf
chr2	11273458	.	A	AC	.	PASS	.
chr2	11273459	.	C	T	.	PASS	.
```

The same records in reverse processing order:

```vcf
chr2	11273459	.	C	T	.	PASS	.
chr2	11273458	.	A	AC	.	PASS	.
```

## Command

```bash
java -Xmx8g -jar snpEff.jar \
  -noLog -noStats -nodownload \
  -c snpEff.config -dataDir /path/to/data \
  GRCh37.75 input.vcf > output.vcf
```

## Observed result

For `chr2:11273459 C>T`:

| input | `WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS` |
|---|---|
| target alone | absent |
| insertion, then target | present |
| target, then insertion | absent |

The token is added to the `C2orf50` consequences for `ENST00000381585` and
`ENST00000405022`. Consequence, impact, gene, transcript, rank and HGVS values
are otherwise identical.

## Expected result

The transcript sanity warning for a variant should not depend on whether an
adjacent variant was processed earlier in the same VCF. The target's `ANN`
value should be consistent across these inputs.

## Environment

- SnpEff 5.4c (2026-02-23)
- Database `GRCh37.75`
- Eclipse Temurin OpenJDK 21.0.12.1
- Linux x86_64, kernel 7.0.0-14-generic
- JAR SHA-256: `5e8f75cbf908a33c6fb2e65c81e66fe31236cb21bb0195541c4703d8202c22b3`

Could a transcript object or a lazily cached protein/CDS value be mutated while
processing the first indel and then reused for the following record? I can
provide the three complete minimal VCFs and resulting outputs if helpful.
````
