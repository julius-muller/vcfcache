#!/usr/bin/env bash
#SBATCH --job-name=vcfcache-pair
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --output=/results/slurm-%A_%a.out

set -euo pipefail

repo_root=${VCFCACHE_REPO_ROOT:?VCFCACHE_REPO_ROOT is required}
task_manifest=${VCFCACHE_TASK_MANIFEST:?VCFCACHE_TASK_MANIFEST is required}
results_root=${VCFCACHE_RESULTS_ROOT:-/results}
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
cache=/mnt/data/vcfcache_data/caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything
params=$cache/params.snapshot.yaml
run_root=/mnt/data/vcfcache_runs/task-${task_id}
result_dir=$results_root/tasks/task-${task_id}
partial_result=${result_dir}.partial-${SLURM_JOB_ID}

if [[ -d "$result_dir" ]]; then
  printf 'Task %s is already archived at %s\n' "$task_id" "$result_dir"
  exit 0
fi

row=$(awk -F '\t' -v row="$((task_id + 2))" 'NR == row { print; exit }' "$task_manifest")
if [[ -z "$row" ]]; then
  printf 'No task row %s in %s\n' "$task_id" "$task_manifest" >&2
  exit 1
fi

IFS=$'\t' read -r manifest_task_id sample input_vcf replicate first_mode second_mode randomization_key <<<"$row"
if [[ "$manifest_task_id" != "$task_id" ]]; then
  printf 'Manifest task mismatch: expected %s, found %s\n' "$task_id" "$manifest_task_id" >&2
  exit 1
fi

rm -rf -- "$run_root"
mkdir -p "$run_root"
test -x "$repo_root/.venv/bin/python"
test -f "$repo_root/benchmarks/run_pilot.py"
export LC_ALL=C
export LANG=C

runner=(
  "$repo_root/.venv/bin/python"
  "$repo_root/benchmarks/run_pilot.py"
)
common=(
  --data-root "$run_root"
  --input "$input_vcf"
  --cache "$cache"
  --params "$params"
  --replicate "$replicate"
)

"${runner[@]}" preflight "${common[@]}"
"${runner[@]}" run --mode "$first_mode" "${common[@]}"
"${runner[@]}" run --mode "$second_mode" "${common[@]}"
"${runner[@]}" compare "${common[@]}"
"${runner[@]}" summarize "${common[@]}"

pilot_dir="$run_root/pilot/$sample"
summary=$(find "$pilot_dir" -name "summary_r$(printf '%02d' "$replicate").json" -print -quit)
comparison=$(find "$pilot_dir" -name "semantic_comparison_r$(printf '%02d' "$replicate").json" -print -quit)
test -n "$summary"
test -n "$comparison"
"$repo_root/.venv/bin/python" - "$comparison" <<'PY'
import json
import sys
from pathlib import Path

report = json.loads(Path(sys.argv[1]).read_text())
if report.get("semantic_pass") is not True:
    raise SystemExit("Semantic comparison did not pass")
PY

rm -rf -- "$partial_result"
mkdir -p "$partial_result"
rsync -a "$run_root/" "$partial_result/"
cat >"$partial_result/slurm-task.json" <<EOF
{
  "task_id": $task_id,
  "sample": "$sample",
  "replicate": $replicate,
  "first_mode": "$first_mode",
  "second_mode": "$second_mode",
  "randomization_key": "$randomization_key",
  "slurm_job_id": "$SLURM_JOB_ID",
  "slurm_node": "$SLURMD_NODENAME"
}
EOF
mv "$partial_result" "$result_dir"
printf 'Archived task %s at %s\n' "$task_id" "$result_dir"
