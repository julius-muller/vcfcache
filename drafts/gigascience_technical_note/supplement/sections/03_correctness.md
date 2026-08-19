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

All 156 VEP and 156 fastVEP cached external-WGS outputs passed. The documented
VEP exception comprised 8,898 HGNC_ID-only record differences among
734,582,949 compared records (0.0012%); the maximum was 495 records in one
4.84-million-variant comparison. A separate 89-output technical validation set
contained 265 HGNC_ID-only differences and no unexpected annotation mismatch.
The 12 cached pipeline-complexity outputs also passed their specified gates.
These results support equivalence within a fixed annotator and recipe, not
between VEP and fastVEP or between arbitrary annotation configurations.
