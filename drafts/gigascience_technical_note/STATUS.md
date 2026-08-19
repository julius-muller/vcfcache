# Draft status

Status: publication-format working draft

## Complete

- GigaScience Technical Note section structure
- frozen claim ledger
- three journal-style main figures with separately exported panels
- real-world and controlled-engineering supplementary figures
- figure provenance, source tables, and checksums
- concise abstract, Findings, Methods, availability, and declarations
- reproducible Markdown, LaTeX, DOCX, and PDF build with justified body text
- independent-sample VEP and fastVEP performance results
- complete 72-cell matched Panel/WES/WGS extension with 48 newly measured
  Panel/WES cells and no execution or semantic failures
- held-out-WGS fastVEP CPU-by-annotation-complexity experiment with strict
  output validation
- explicit correctness accounting, including the 265 HGNC_ID-only differences

## Next editorial pass

1. Confirm title, author list, affiliations, and corresponding author.
2. Complete scientific and editorial review of the selected three main and three
   supplementary figures.
3. Complete the code/data archive DOI and software identifiers.
4. Obtain the institutionally appropriate public-data ethics wording.
5. Resolve all items in `OPEN_ITEMS.md`, then run `make release-check`.

Alternative title retained for discussion: *Annotate what is new: VCFcache
accelerates repeated variant annotation without changing results*.
