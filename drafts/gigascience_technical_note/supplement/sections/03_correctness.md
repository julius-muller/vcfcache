## Semantic validation details

Correctness was a gate, not a secondary timing outcome. A cached result could
not contribute to a performance figure unless it passed the comparator defined
for its annotator.

For VEP, variant keys, retained input INFO, genotype, CSQ header and all CSQ
entries were compared. CSQ-entry order and complete-record order within a
shared locus were canonicalized. HGNC_ID-only differences were counted but
excluded from pass/fail because VEP 115.2 exhibits a documented batching
dependency. This exception was not generalized to other fields or annotators.

For fastVEP, every INFO and FORMAT value and every relevant header definition
was compared. INFO-tag order, CSQ-entry order and complete-record order within
identical loci were canonicalized; no annotation field was ignored. This
included fastVEP supplementary fields such as FV_CLINVAR, FV_GNOMAD, FV_1KG
and FV_TOPMED when present.

For SnpEff, the comparator checked variant keys, genotype and the complete INFO
content defining the tested recipe. ANN entries and same-position records were
canonicalized for ordering; warning tokens, ClinVar and gnomAD values were not
ignored. Of 180 final-campaign pairs, 176 passed. Four cells from GRCh37 sample
`hu6664F3` were excluded before inference because one allele carried a
`WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS` token only in the cached arm. All
consequence, transcript, HGVS, ClinVar, gnomAD, genotype and record-count values
otherwise agreed.

A direct-SnpEff-only diagnostic reduced the difference to two input records.
For `chr2:11273459 C>T`, the warning was absent when the allele was processed
alone, present after `chr2:11273458 A>AC`, and absent when the same two records
were reversed. The target's warning-free ANN digest was identical in both
controls. This demonstrates processing-order-dependent SnpEff 5.4c diagnostic
state rather than a VCFcache lookup or merge error. The warning was neither
suppressed nor added to a general exception; exact VCFs and output checksums are
included in the archived diagnostic.

All 156 VEP and 156 fastVEP cached external-WGS outputs passed. The documented
VEP exception comprised 8,898 HGNC_ID-only record differences among
734,582,949 compared records (0.0012%); the maximum was 495 records in one
4.84-million-variant comparison. A separate 89-output technical validation set
contained 265 HGNC_ID-only differences and no unexpected annotation mismatch.
The 12 cached pipeline-complexity outputs also passed their specified gates.
These results support equivalence within a fixed annotator and recipe, not
between annotators or arbitrary annotation configurations.
