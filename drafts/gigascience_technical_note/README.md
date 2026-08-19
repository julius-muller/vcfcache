# VCFcache GigaScience Technical Note

This directory is the independent source workspace for a concise GigaScience
Technical Note and its bioRxiv preprint. It is intentionally separate from the
historical manuscript draft and from benchmark rendering directories.

The Markdown files in `sections/` are the editable source of truth. Figures
are copied from frozen benchmark snapshots and accompanied by checksums and
source paths. Numeric claims are registered in `claims.tsv` before they enter
the prose.

## Working order

1. Edit one section at a time under `sections/`.
2. Record unresolved publication metadata in `OPEN_ITEMS.md`.
3. Run `make check` for structural and provenance checks.
4. Run `make draft` to create combined Markdown, editable DOCX, LaTeX, and PDF.
5. Run `make release-check` only when all author-controlled placeholders have
   been resolved.

The draft targets a 180–220 word structured abstract, approximately 1,200–
1,600 words of Findings, and approximately 900–1,400 words of Methods. The
supplement holds detailed cohort preparation, execution, and validation
protocols that are needed for reproducibility but not for the central story.

## Evidence boundary

Primary performance claims use only the complete 52-genome VEP/fastVEP
snapshot dated 2026-08-09 and the final dual-cache pipeline-complexity snapshot
dated 2026-08-12. KPGP, SGDP and PGP evaluation genomes were not used to build
the cohort caches and had no documented project overlap with gnomAD.

Main Figure 2 and Supplementary Figure 3 add a matched assay-scale extension
using 12 independent GRCh38 KPGP and SGDP genomes, the bundled gnomAD AF ≥1%
strategy, and exactly one retained direct and cached observation per
sample-tool-assay cell. No 1000 Genomes-derived performance sample is used.

Supplementary Figure 2 adds a controlled computing interaction experiment on
one held-out KPGP genome. It compares core and deliberately dense fastVEP
recipes under process-wide CPU affinity at a fixed approximately 90% hit rate;
it is not used to estimate biological cache coverage.

The older assay-size campaign based on 1000 Genomes-derived inputs is excluded
from performance figures. Its semantic-comparison audit is retained only as a
technical correctness record.

## Outputs

Generated files are placed in `build/` and are not the editable source:

- `vcfcache_gigascience_technical_note.md`
- `vcfcache_gigascience_technical_note.docx`
- `vcfcache_gigascience_technical_note.tex`
- `vcfcache_gigascience_technical_note.pdf`
- `vcfcache_supplementary_material.md`
- `vcfcache_supplementary_material.docx`
- `vcfcache_supplementary_material.tex`
- `vcfcache_supplementary_material.pdf`

Body paragraphs are justified in the editable DOCX and review PDF. Equations
are authored as LaTeX mathematics in Markdown and are emitted as native Word
equations and as journal-portable LaTeX source.

A future journal-specific draft should be created as a sibling directory, not
by changing this workspace into a different article type.
