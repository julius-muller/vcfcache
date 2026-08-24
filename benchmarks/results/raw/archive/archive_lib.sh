#!/bin/bash
# Shared helpers so every campaign lands in the archive by construction.
#
# Runs are never deleted or overwritten: each writes into its own directory
# keyed by configuration, and a repeat of the same configuration is suffixed
# with its Slurm job id rather than replacing what is already there.
ARCHIVE_ROOT=${ARCHIVE_ROOT:-/results/benchmark_archive}

archive_dir () {  # archive_dir <run_id> -> prints the directory to write into
  local id="$1" dir="$ARCHIVE_ROOT/$1"
  if [[ -s "$dir/metrics.tsv" || -d "$dir/tasks" ]]; then
    dir="$ARCHIVE_ROOT/${id}-job${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
  fi
  mkdir -p "$dir"
  printf '%s' "$dir"
}

write_manifest () {  # write_manifest <dir> <json>
  # Only the first array task writes it; the rest would race on the same file.
  local dir="$1" json="$2"
  [[ "${SLURM_ARRAY_TASK_ID:-0}" == "0" ]] || return 0
  printf '%s\n' "$json" > "$dir/RUN_MANIFEST.json"
}
