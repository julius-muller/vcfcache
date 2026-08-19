# Data availability

The datasets supporting the benchmark results are available from
https://doi.org/10.5281/zenodo.22018995. The deposit will include frozen
source tables, cohort-selection and QC manifests, commands, semantic-comparison
summaries, figure-rendering code, R session records, and a compact test input
with expected cached and direct outputs.

A self-contained executable example is distributed with the software itself and
requires no download, no reference data and no external annotation tool. The
command `vcfcache demo` runs the complete workflow — blueprint creation,
blueprint extension, cache construction and cached annotation — against bundled
test data using a self-describing annotation recipe, then asserts that the
cached and directly annotated outputs are identical. It completes in
approximately one second and its peak resident memory is 42.6 MB, so editors,
reviewers and users can verify the central correctness claim on any machine
that satisfies the stated requirements. This example is the recommended entry
point for reproducing the workflow; the archived benchmark deposit reproduces
the measurements. Raw public human genomic inputs will
be referenced through their provider accessions and terms rather than
redistributed when redistribution is not appropriate.

The bundled GRCh38 VEP cache with any-stratum AF ≥10% is archived at
https://doi.org/10.5281/zenodo.18189447, and the AF ≥1% cache at
https://doi.org/10.5281/zenodo.18190046. The corresponding GRCh37 caches are
archived at https://doi.org/10.5281/zenodo.18189051 and
https://doi.org/10.5281/zenodo.18189348. Source software will be cited by
release and persistent archive identifier in the final manuscript.
