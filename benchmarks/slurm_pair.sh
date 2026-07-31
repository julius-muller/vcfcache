#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --exclusive
#SBATCH --time=12:00:00

set -euo pipefail

repo_root=${VCFCACHE_REPO_ROOT:?VCFCACHE_REPO_ROOT is required}
task_manifest=${VCFCACHE_TASK_MANIFEST:?VCFCACHE_TASK_MANIFEST is required}
results_root=${VCFCACHE_RESULTS_ROOT:-/results}
campaign_id=${VCFCACHE_CAMPAIGN_ID:?VCFCACHE_CAMPAIGN_ID is required}
phase=${VCFCACHE_PHASE:?VCFCACHE_PHASE is required}
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
cache=/mnt/data/vcfcache_data/caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything
params=$cache/params.snapshot.yaml
run_root=/mnt/data/tmp/vcfcache_runs/$campaign_id/$phase/task-$task_id
phase_root=$results_root/campaigns/$campaign_id/$phase
result_dir=$phase_root/tasks/task-$task_id
attempt_dir=$phase_root/attempts/task-$task_id/job-$SLURM_JOB_ID
partial_result=$result_dir.partial-$SLURM_JOB_ID
completed=0

archive_attempt() {
  local status=$?
  if ((completed)); then
    return 0
  fi
  mkdir -p "$attempt_dir"
  if [[ -d "$run_root" ]]; then
    rsync -a "$run_root/" "$attempt_dir/run/"
  fi
  cat >"$attempt_dir/failure.json" <<EOF
{
  "campaign_id": "$campaign_id",
  "phase": "$phase",
  "task_id": $task_id,
  "slurm_job_id": "$SLURM_JOB_ID",
  "slurm_array_job_id": "${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}",
  "slurm_node": "${SLURMD_NODENAME:-unknown}",
  "exit_code": $status
}
EOF
  return "$status"
}
trap archive_attempt EXIT

if [[ -d "$result_dir" ]]; then
  printf 'Task %s/%s/%s is already archived at %s\n' \
    "$campaign_id" "$phase" "$task_id" "$result_dir"
  completed=1
  exit 0
fi

case "$campaign_id:$phase:$task_id" in
  *[!A-Za-z0-9_.:-]*)
    printf 'Unsafe campaign, phase, or task identifier.\n' >&2
    exit 1
    ;;
esac

row=$(awk -F '\t' -v row="$((task_id + 2))" 'NR == row { print; exit }' "$task_manifest")
if [[ -z "$row" ]]; then
  printf 'No task row %s in %s\n' "$task_id" "$task_manifest" >&2
  exit 1
fi

IFS=$'\t' read -r \
  manifest_task_id manifest_phase measured sample population superpopulation sex \
  input_vcf input_records input_sha256 replicate first_mode second_mode \
  randomization_key <<<"$row"
if [[ "$manifest_task_id" != "$task_id" || "$manifest_phase" != "$phase" ]]; then
  printf 'Manifest mismatch for task %s/%s.\n' "$phase" "$task_id" >&2
  exit 1
fi
if [[ "$first_mode:$second_mode" != "cached:uncached" && \
      "$first_mode:$second_mode" != "uncached:cached" ]]; then
  printf 'Invalid paired mode order in task manifest.\n' >&2
  exit 1
fi

rm -rf -- "$run_root"
mkdir -p "$run_root" "$attempt_dir"
test -x "$repo_root/.venv/bin/python"
test -f "$repo_root/benchmarks/run_pilot.py"
test -f "$input_vcf"
test -f "$input_vcf.tbi"
export LC_ALL=C
export LANG=C
export VCFCACHE_REQUIRE_CGROUP_METRICS=1

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
srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=8 \
  "${runner[@]}" run --mode "$first_mode" "${common[@]}"
srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=8 \
  "${runner[@]}" run --mode "$second_mode" "${common[@]}"
"${runner[@]}" compare "${common[@]}"
"${runner[@]}" summarize "${common[@]}"

pilot_dir=$run_root/pilot/$sample
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

mkdir -p "$phase_root/tasks"
test ! -e "$partial_result"
mkdir "$partial_result"
rsync -a "$run_root/" "$partial_result/"
cat >"$partial_result/slurm-task.json" <<EOF
{
  "campaign_id": "$campaign_id",
  "phase": "$phase",
  "measured": $measured,
  "task_id": $task_id,
  "sample": "$sample",
  "population": "$population",
  "superpopulation": "$superpopulation",
  "sex": "$sex",
  "input_vcf": "$input_vcf",
  "input_records": $input_records,
  "input_sha256": "$input_sha256",
  "replicate": $replicate,
  "first_mode": "$first_mode",
  "second_mode": "$second_mode",
  "randomization_key": "$randomization_key",
  "slurm_job_id": "$SLURM_JOB_ID",
  "slurm_array_job_id": "${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}",
  "slurm_node": "$SLURMD_NODENAME",
  "allocated_cpus": "${SLURM_CPUS_PER_TASK:-8}",
  "requested_memory_per_node": "${SLURM_MEM_PER_NODE:-30000}"
}
EOF
mv "$partial_result" "$result_dir"
cat >"$attempt_dir/success.json" <<EOF
{
  "campaign_id": "$campaign_id",
  "phase": "$phase",
  "task_id": $task_id,
  "slurm_job_id": "$SLURM_JOB_ID",
  "result_dir": "$result_dir",
  "exit_code": 0
}
EOF
completed=1
printf 'Archived task %s/%s/%s at %s\n' \
  "$campaign_id" "$phase" "$task_id" "$result_dir"
