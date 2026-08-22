# Benchmarking guidelines

Rules that earned their place. Most were written after a run had to be thrown
away or repeated.

## Store the file, not the number

If a tool writes something, keep what it wrote. A pipeline that ends
`| awk '{print $3}' >> summary.tsv` has destroyed the evidence and kept an
opinion about it.

```
# wrong: only the number survives, and only this question can ever be asked
bcftools index -s in.bcf | awk '{s+=$3} END {print s}' > n_records.txt

# right: keep the tool's output, derive the number when you need it
bcftools index -s in.bcf > raw/in.bcf.index_stats.txt
```

The same applies to `bcftools stats`, `/usr/bin/time -v`, annotator stderr,
Slurm `.out` files and cache-build logs. Storage is cheap; a rerun on a shared
cluster is not, and a reviewer's question you cannot answer is worse.

Derived tables are fine and expected — they just never replace their source.
Every derived value should be reproducible from something in `raw/`.

## Measure the thing you are claiming

- **Separate start-up from per-variant work.** Caching cannot touch start-up, so
  a single wall time cannot tell you whether caching helped for the right
  reason. Measure at two input sizes and fit `t = s + c*N`.
- **Never size from one small slice.** An 11-hour job was once planned from a
  ratio read as a rate. Two points, or no estimate.
- **Both arms in the same execution shape.** Splitting an annotation by contig
  gives the direct arm a parallel speedup the smaller cached arm cannot use.
  Measured on the same genomes: serial 6.85x, split 3.09x. Timed arms are always
  the documented single-process invocation.
- **Warm-up belongs to short runs only.** Cold state dominates a sub-minute
  Panel or WES cell, so those get one untimed pre-run. At WGS scale it costs an
  entire extra pass and changes nothing measurable.
- **Condition order is a confounder.** One Panel/WES campaign had to be excluded
  entirely because cold-state order, not caching, explained the short cells.

## Assert configuration, never infer it

Every silent failure so far came from an unasserted assumption.

- A hardcoded genome name that no longer existed after a clean reinstall.
- A record count parsed from a command that prints a warning to **stdout**, so
  the count was non-numeric, evaluated as zero, and silently took the wrong
  branch. Use `bcftools index -s | awk '{s+=$3}'` and assert the result is a
  number.
- A `params.yaml` built from the GRCh38 template for a GRCh37 run, with only
  some fields rewritten.
- A pipeline stage that died mid-stream and left a well-formed but
  under-annotated BCF, which would have been timed as if it were complete.

So: every runner logs the path it took, fails loudly on a bad count, and asserts
that the fields defining the level are actually present in its output.

## Gate before timing

A new cache is proven on a chr22 slice — direct versus cached, comparing every
field the level is defined by — before any measurement is recorded against it.
A cache built from a wrong recipe produces perfectly consistent, perfectly
meaningless numbers.

## Comparing annotation output

Records sharing a `CHROM`+`POS` can come back in a different order from a
positional annotator (multi-allelic split order, indel-versus-SNV tie-break)
without any annotation differing. A plain order-sensitive digest therefore fails
on essentially every real sample.

Do not reach for `sort | md5sum`: over 4.6 M records it took longer than the
annotation itself and ran every task into the wall clock limit. Sort **within**
each position block instead — an annotator never moves a record across
positions — which is a single streaming pass. See `raw/` runs from 2026-08-21
onward; the block-local digest cost 29 s where the external sort exceeded hours.

Beware `set -o pipefail` with `grep -q`: `grep` closes the pipe on its first
match, the producer dies of SIGPIPE, and the pipeline reports failure *because*
the pattern was found. This killed a 40-minute run that had worked.

## Never discard a successful run

Every run that completed goes into the archive under an id describing its
configuration — `<tool>-<level>-<assembly>-<assays>-<invocation>` — with a
`RUN_MANIFEST.json` recording tool versions, database, cache, Slurm job id and
status (`primary`, `superseded`, `diagnostic`). A repeat of the same
configuration is suffixed with its job id rather than overwriting.

Superseded runs are kept, not deleted: the split-invocation SnpEff runs are the
optimised-invocation lower bound and answer a question the serial runs cannot.

## Getting numbers into the manuscript

A number reaches the text only through
`drafts/gigascience_technical_note/claims.tsv`, with its evidence file present,
and with `make check` passing. If the evidence file is missing, the claim is not
ready.
