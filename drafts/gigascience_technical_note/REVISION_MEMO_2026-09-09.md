# Response to the Elicit citation and revision memo

Source: author-supplied `Elicit - VCFcache technical-note citation and revision memo.md`.
Reviewed against the current section sources and `vcfcache/database/workflow_manager.py`.

The central story remains reusable computation across samples, with the measured
9.50x VEP, 2.54x clinical fastVEP and 10.78x clinical SnpEff WGS gains,
cohort build amortization, and absolute time returned. The existing title is
retained. Precision about the execution contract strengthens these claims.

| Feedback | Revision or disposition |
|---|---|
| Scope novelty to incremental execution | Findings now identifies payload reuse, automatic miss execution and sample output assembly as the contribution. Adjacent lookup and cohort systems are discussed without claiming a comprehensive novelty search. |
| Avoid “complete per-allele output” | Abstract and Findings use INFO annotation payloads; Methods defines preserved sample data and unsupported arbitrary record/FORMAT transformations. |
| Correct perfect-hit runtime | Separate direct/cached non-annotation costs and conditional startup; specify the equal-average-cost assumption for hits and misses and explain perfect-hit break-even. |
| Interpret ideal figure curves | Captions describe the curves as annotation-dominated approximations. Existing plotted observations and exports are preserved. |
| Specify normalization/key | State CHROM/POS/REF/ALT, assembly scope, split records, ALT=* exclusion, contig compatibility and upstream reference normalization. No claim that splitting alone left-aligns. |
| Specify recipe fingerprint | Distinguish product snapshots/YAML identity from benchmark resource checksums; identify recipe changes requiring a new cache. |
| Specify equivalence | Keep all existing comparator definitions and failures; explicitly state that the HGNC exception and excluded SnpEff cells do not establish equality of those outputs. |
| Separate population data from cache provenance | Cite gnomAD, its release and the existing Zenodo cache artifacts separately; explain cohort-specific coverage and local fastVEP/SnpEff annotation. |
| Add stage profile | Supplement reports a real-WGS command-level breakdown at 92.99% hits recovered from the archived measured run, alongside the separate controlled perfect-hit diagnostics. No new benchmark is required. |
| Broaden citations | Add Tan, GEMINI, VarNote and gnomAD. VRS and GenomicsDB are omitted to keep the note concise; VCFcache does not implement VRS identifiers. |

## Citation corrections

The memo's VarNote DOI was incorrect: the article uses
[10.1101/gr.267997.120](https://genome.cshlp.org/content/30/12/1789),
not `10.1101/gr.264192.120`. GEMINI's third author is **Rory Kirchner**,
not Ryan. These were checked against the article records. Tan's bibliographic
record is [10.1093/bioinformatics/btv112](https://pubmed.ncbi.nlm.nih.gov/25701572/);
the gnomAD resource is [10.1038/s41586-020-2308-7](https://www.nature.com/articles/s41586-020-2308-7).
The existing peer-reviewed echtvar reference is retained rather than substituting
the memo's preprint footnote.

## Remaining evidence work

No new benchmark execution is implied by this editorial pass. Existing measured
headline claims remain unchanged. The stage breakdown is now documented from
existing evidence. The normalization provenance audit remains in `OPEN_ITEMS.md`;
author, ethics and archive
metadata remain outstanding. Do not describe the draft as submission-final.

## Correction following evidence recheck

The initial revision incorrectly stated that stage timing with cache misses
remained to be measured. Both the early pilot's `stage_timings.tsv` (80%, 90%
and 100% conditions) and real-WGS `workflow.log` files already record command
durations. The supplement now includes a measured real-WGS example reconciled
to its retained metrics row. None of the 180 cached final-campaign cells has
zero misses. The 100% diagnostic belongs to a constructed engineering cache.
Normalization wording was also refined: consistent allele keys enable reuse;
left alignment is not an automatic consequence of multiallelic splitting or
a prerequisite for matching records that already share the same representation.
