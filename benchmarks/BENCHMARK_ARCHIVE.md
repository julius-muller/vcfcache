# Benchmark archive

Every successful benchmark run is stored permanently under
`/results/benchmark_archive` on the cluster, one directory per run. Nothing is
deleted or overwritten, including runs superseded by a later configuration: a
superseded measurement still bounds what that configuration achieves, and is
needed to show, for example, what an optimised SnpEff invocation reaches
relative to the documented serial one.

## Layout

    <run_id>/
      metrics.tsv        one row per measured cell
      RUN_MANIFEST.json  tool, version, database, recipe, invocation, assembly,
                         annotation level, slurm job, and how it was produced
      provenance.txt     source paths the metrics were derived from

`archive_lib.sh` in the same directory provides `archive_dir` and
`write_manifest`. Campaigns source it and write into the archive by
construction; a repeat of an existing configuration is suffixed with its Slurm
job id rather than replacing what is there.

## Naming

    <tool>-<level>-<assembly>-<assays>-<invocation>

`level` is `core` or `clinical` (see `ANNOTATION_LEVELS.md`). `invocation`
records how the annotator was driven, because it materially changes the result:
`serial` is the documented single-process form, `split` runs one annotator
process per contig.

## Status

`RUN_MANIFEST.json` carries `status`:

- `primary` — intended for the manuscript
- `superseded` — a correct measurement of a configuration no longer primary
- `diagnostic` — calibration or probe, not a campaign result

## Recovering a lost run

Job logs under `/results/<campaign>/logs/` are retained and carry the timed
lines, so a metrics file deleted by mistake can be reconstructed. The
`snpeff-core-grch38-wgs-split` run was rebuilt this way; its manifest records
that it is a reconstruction rather than an original file.
