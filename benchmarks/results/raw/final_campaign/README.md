# Canonical final-campaign index

This directory is the canonical index for the campaign reported in
`../../../RESULTS.md`. The untouched task TSVs remain in `../archive/`.

Slurm array tasks started concurrently. The original `archive_dir` helper saw
the first task's new directory and routed the other eleven tasks to a
`-job<id>` directory. No task was lost: every arm has one task TSV in the base
directory and eleven in the job-suffixed directory. `RUN_MANIFEST.json` lists
both source directories for each arm and corrects the array-task manifest's
misleading single-cohort field without altering the original evidence.

The accepted runner is frozen under
`../../../reproduction/final_clinical_campaign/`. Its GRCh37 input discovery
used `ls | head -12`; `samples.tsv` freezes the exact selected files and must be
used instead in any future run.
