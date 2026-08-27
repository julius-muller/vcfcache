# hu6664F3 SnpEff concordance diagnostic

Untimed retained-output job 1205 repeated the GRCh37 WES cell with the exact
core and clinical recipes. Direct and cached outputs contained the same 86,997
variant keys. Each level had exactly one differing record:
`chr2:11273459 C>T`.

For two `C2orf50` 5-prime UTR transcript consequences, the cache-build SnpEff
annotation includes `WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS`, whereas direct
whole-exome SnpEff annotation omits it. Consequence, impact, gene, transcript,
HGVS, ClinVar and gnomAD AF values otherwise agree. The clinical and core
comparisons identify the same warning-only difference; list-order
canonicalization does not remove it.

Direct-SnpEff-only jobs 1236 and 1238 subsequently established a minimal
processing-order reproducer. Annotating `chr2:11273459 C>T` alone emits no
warning. Annotating the adjacent `chr2:11273458 A>AC` first and then the target
adds `WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS` to the target's two `C2orf50`
consequences. Reversing those same two records so that the target is processed
first removes the warning again. The target's remaining `ANN` content is
unchanged.

This is direct SnpEff 5.4c behavior with the same `GRCh37.75` database and no
VCFcache process in the path. It proves an unintended SnpEff processing-order
side effect rather than use of another variant's INFO values as annotation
evidence or a VCFcache merge error. The four original WES/WGS cells remain
excluded from strict inferential summaries because the raw records are not
identical. No broad SnpEff warning exception is added to the publication
comparator.

The complete BCFs remain on the Slurm result volume at
`/mnt/data/slurm-results/diagnostics/hu6664F3-1205/`. The TSVs here contain the
record counts and exact differing annotation values; `*.sha256.txt` freezes the
retained BCF identities.

The complete context matrix, minimal two-record inputs and direct SnpEff outputs
are frozen under `snpeff_context/`. `context_matrix.tsv` shows that the sample
context remains warning-free from one record through the complete chromosome,
whereas the blueprint changes from warning-free at one record to warning-bearing
within the two-kilobase window. `minimal_matrix.tsv` isolates processing order.
