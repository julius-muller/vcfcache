# Final native-annotation figure set

This snapshot is additive and does not overwrite any prior figure set.
All plots were rendered in R with ggplot2; composite panels use patchwork.
PNG files are 600 dpi, and every manuscript panel is also exported as SVG and PDF.

## Manuscript figures

- `figure_assay_three_annotators`: matched Panel/WES/WGS results for 12 held-out GRCh38 genomes. VEP `--everything`, clinical fastVEP and clinical SnpEff are each compared only with their own direct baseline; their annotation outputs are not claimed to be identical across tools.
- `figure_snpeff_annotation_load`: within-SnpEff GRCh37 core-versus-clinical contrast on the same samples and cache membership.
- `figure_pipeline_time_returned`: controlled VEP annotation-load experiment at 80.23% and 90.26% observed cache coverage.

## Repository and alternative figures

- `headliner_native_wgs`: deliberately simple normalized-time view for the repository front page.
- `native_wgs_clock_lines`: compact direct-to-cached runtime alternative.

## Interpretation boundary

The combined annotator panels show within-tool speedups under native recipes. They must not be used to rank annotation completeness or claim cross-tool semantic equivalence. The 12 genomes are independent KPGP/SGDP evaluation samples and are not 1000 Genomes samples.

