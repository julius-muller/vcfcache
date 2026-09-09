# Draft status

Status: publication-format working draft

## Complete

- September 9 citation/implementation review: scoped INFO-payload contribution,
  conditional-startup runtime model, explicit normalization/resource contracts
  and qualified semantic equivalence; see `REVISION_MEMO_2026-09-09.md`

- GigaScience Technical Note section structure
- frozen claim ledger
- four journal-style main figures with separately exported panels
- real-world and controlled-engineering supplementary figures
- figure provenance, source tables, and checksums
- concise abstract, Findings, Methods, availability, and declarations
- reproducible Markdown, LaTeX, DOCX, and PDF build with justified body text
- independent-sample VEP and fastVEP performance results
- final 180-cell fastVEP/SnpEff native annotation campaign, with 176 strict
  passes and four transparently excluded upstream warning-only SnpEff pairs
- direct-SnpEff minimal reproducer and copy-ready upstream issue report
- complete 108-pair matched Panel/WES/WGS comparison spanning VEP, clinical
  fastVEP and clinical SnpEff, with no 1000 Genomes-derived performance sample
- held-out-WGS fastVEP CPU-by-annotation-complexity experiment with strict
  output validation
- explicit correctness accounting, including the 265 HGNC_ID-only differences

## Next editorial pass

1. Confirm title, author list, affiliations, and corresponding author.
2. Complete scientific and editorial review of the selected four main and three
   supplementary figures.
3. Complete the code/data archive DOI and software identifiers.
4. Obtain the institutionally appropriate public-data ethics wording.
5. Resolve all items in `OPEN_ITEMS.md`, then run `make release-check`.
6. Complete the normalization provenance audit. The representative WGS stage
   breakdown has been recovered from archived logs and included in the supplement.

Alternative title retained for discussion: *Annotate what is new: VCFcache
accelerates repeated variant annotation without changing results*.
