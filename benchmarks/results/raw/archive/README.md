# Benchmark archive

Every successful benchmark run is stored here permanently, one directory per
run, and is never deleted or overwritten. A run that was superseded by a later
configuration is still kept: superseded measurements bound what a configuration
can achieve, and are needed to show, for example, what an optimised SnpEff
invocation reaches relative to the documented serial one.

## Layout

    <run_id>/
      metrics.tsv        one row per measured cell
      RUN_MANIFEST.json  tool, version, database, recipe, invocation, assembly,
                         annotation level, slurm job, and how it was produced
      provenance.txt     source paths the metrics were derived from

## Naming

    <tool>-<level>-<assembly>-<assays>-<invocation>

`level` is `core` (consequences only) or `clinical` (core plus ClinVar and
gnomAD allele frequencies). `invocation` records how the annotator was driven,
because that materially changes the result: `serial` is the documented
single-process form; `split` runs one annotator process per contig.

## Status of a run

`RUN_MANIFEST.json` carries `status`:

- `primary`   — intended for the manuscript
- `superseded`— correct measurement of a configuration no longer used as primary
- `diagnostic`— calibration or probe, not a campaign result
