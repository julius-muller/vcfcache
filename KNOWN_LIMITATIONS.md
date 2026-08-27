# Known limitations

Differences between `vcfcache annotate` and running the annotator directly, found
while benchmarking. **None of them change annotation content** — the same variants
get the same annotations — but each one makes a naive `diff` of the two outputs
report a mismatch, so users comparing cached against direct output need to know.

Each entry is a candidate fix, not settled behaviour.

## 1. Contig order follows the cache, not the input

Cached output emits contigs in the **cache's** order; direct annotation preserves
the **input VCF's** order. When the two differ, every record after the first
displaced contig appears at a different line, and a byte-wise comparison fails
even though every record is identical.

Observed on GRCh37 PGP genomes against a GRCh38-ordered cache: 185 aligned keys,
`CSQ`, `FV_CLINVAR` and `FV_GNOMAD` all identical, whole contig blocks displaced.
GRCh38 cohorts were unaffected only because their contig order happened to match.

*Impact*: `diff`/`md5sum` comparisons fail; downstream tools that assume input
contig order may need a re-sort. *Possible fix*: emit in the input's contig order,
or document that output follows the cache and offer a flag.

## 2. Exact duplicate records are collapsed

If an input VCF lists the same `CHROM/POS/REF/ALT` more than once, cached output
keeps one copy; direct annotation keeps them all.

Measured on `hu365511` (PGP GRCh37), whole genome:

```
direct   5,005,430 records   4,985,467 distinct keys
cached   4,989,016 records   4,985,467 distinct keys
```

The distinct sets are identical — only duplicate copies differ.

*Impact*: record counts shrink silently; a user checking `bcftools index -n`
before and after sees fewer variants and may reasonably think data was lost.
*Possible fix*: preserve duplicates, or log a count of collapsed records the way
dropped variants are already logged.

## 3. Variants the annotator drops are removed from cached output

Already deliberate and logged (`workflow_manager.py`, step 4/4: "Annotation tool
dropped N variants from input. By default, these are also removed from cached
output to match annotation tool behavior"). Recorded here so the three behaviours
that change record counts or order sit in one place.

## Comparing cached and direct output

A digest that tolerates both differences above, in a single streaming pass:

```bash
# per-contig digest, combined in sorted contig order; within each CHROM+POS block
# records are sorted, since a positional annotator may reorder those freely
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' out.bcf \
  | awk -f benchmarks/results/raw/blocksort.awk \
  | awk -v d="$tmp" '{ print >> (d"/"$1) }'
for c in $(ls "$tmp" | sort); do printf '%s %s\n' "$c" "$(md5sum < "$tmp/$c" | cut -d' ' -f1)"; done | md5sum
```

Do not use `sort | md5sum` over the whole stream: at WGS scale it takes longer
than the annotation it is checking.
