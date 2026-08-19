#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --exclusive
#SBATCH --time=04:00:00

set -euo pipefail

repo_root=${VCFCACHE_REPO_ROOT:?required}
runner=${VCFCACHE_ASSAY_RUNNER:?required}
task_manifest=${VCFCACHE_ASSAY_TASK_MANIFEST:?required}
campaign_id=${VCFCACHE_CAMPAIGN_ID:?required}
task_id=${SLURM_ARRAY_TASK_ID:?required}
results_root=${VCFCACHE_RESULTS_ROOT:-/results}
run_root=/mnt/data/tmp/vcfcache_external_assay/$campaign_id/task-$task_id
result_dir=$results_root/campaigns/$campaign_id/tasks/task-$task_id
attempt_dir=$results_root/campaigns/$campaign_id/attempts/task-$task_id/job-$SLURM_JOB_ID
partial_result=$result_dir.partial-$SLURM_JOB_ID
completed=0

archive_attempt() {
  local status=$?
  if ((completed)); then return 0; fi
  mkdir -p "$attempt_dir"
  [[ ! -d "$run_root" ]] || rsync -a "$run_root/" "$attempt_dir/run/"
  printf '{"exit_code": %d, "job_id": "%s", "task_id": %d}\n' \
    "$status" "$SLURM_JOB_ID" "$task_id" >"$attempt_dir/failure.json"
  return "$status"
}
trap archive_attempt EXIT

if [[ -s "$result_dir/external_assay_summary.json" ]]; then
  completed=1
  exit 0
fi
rm -rf -- "$run_root"
mkdir -p "$run_root" "$(dirname "$result_dir")" "$attempt_dir"
export LC_ALL=C
export LANG=C
export VCFCACHE_REQUIRE_CGROUP_METRICS=1
export PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}"

base=(
  "$repo_root/.venv/bin/python" "$runner"
  --task-manifest "$task_manifest" --task-id "$task_id" --run-root "$run_root"
)
srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=8 \
  "${base[@]}" --condition warmup
mapfile -t conditions < <("${base[@]}" --print-order)
[[ ${#conditions[@]} -eq 2 ]]
for condition in "${conditions[@]}"; do
  srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=8 \
    "${base[@]}" --condition "$condition"
done
srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=3 \
  "${base[@]}" --condition finalize

test -s "$run_root/external_assay_summary.json"
mkdir "$partial_result"
rsync -a "$run_root/" "$partial_result/"
mv "$partial_result" "$result_dir"
printf '{"exit_code": 0, "job_id": "%s", "task_id": %d}\n' \
  "$SLURM_JOB_ID" "$task_id" >"$attempt_dir/success.json"
completed=1
