#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --exclusive
#SBATCH --time=12:00:00

set -euo pipefail

repo_root=${VCFCACHE_REPO_ROOT:?VCFCACHE_REPO_ROOT is required}
task_manifest=${VCFCACHE_TASK_MANIFEST:?VCFCACHE_TASK_MANIFEST is required}
strategies=${VCFCACHE_STRATEGIES:?VCFCACHE_STRATEGIES is required}
results_root=${VCFCACHE_RESULTS_ROOT:-/results}
campaign_id=${VCFCACHE_CAMPAIGN_ID:?VCFCACHE_CAMPAIGN_ID is required}
phase=${VCFCACHE_PHASE:?VCFCACHE_PHASE is required}
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
run_root=/mnt/data/tmp/vcfcache_external_runs/$campaign_id/$phase/task-$task_id
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
  python3 - "$attempt_dir/failure.json" "$status" <<'PY'
import json
import os
import sys
from pathlib import Path

Path(sys.argv[1]).write_text(json.dumps({
    "campaign_id": os.environ["VCFCACHE_CAMPAIGN_ID"],
    "phase": os.environ["VCFCACHE_PHASE"],
    "task_id": int(os.environ["SLURM_ARRAY_TASK_ID"]),
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "slurm_node": os.environ.get("SLURMD_NODENAME", "unknown"),
    "exit_code": int(sys.argv[2]),
}, indent=2, sort_keys=True) + "\n")
PY
  return "$status"
}
trap archive_attempt EXIT

if [[ -d "$result_dir" ]]; then
  completed=1
  exit 0
fi
case "$campaign_id:$phase:$task_id" in
  *[!A-Za-z0-9_.:-]*) exit 1 ;;
esac

rm -rf -- "$run_root"
mkdir -p "$run_root" "$attempt_dir"
export LC_ALL=C
export LANG=C
export VCFCACHE_REQUIRE_CGROUP_METRICS=1

runner=(
  "$repo_root/.venv/bin/python" "$repo_root/benchmarks/run_external_task.py"
  --task-manifest "$task_manifest" \
  --strategies "$strategies" \
  --task-id "$task_id" \
  --run-root "$run_root"
)
mapfile -t conditions < <("${runner[@]}" --print-order)
if ((${#conditions[@]} != 4)); then
  echo "Expected four publication conditions, found ${#conditions[@]}" >&2
  exit 1
fi
for condition in "${conditions[@]}"; do
  srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=8 \
    "${runner[@]}" --condition "$condition"
done

# Comparisons and hashing are deliberately outside all timed annotation cells.
srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=3 \
  "${runner[@]}" --finalize

test -s "$run_root/external_summary.json"
test -d "$phase_root/tasks"
test -d "$phase_root/attempts"
test ! -e "$partial_result"
mkdir "$partial_result"
rsync -a "$run_root/" "$partial_result/"
python3 - "$partial_result/slurm-task.json" <<'PY'
import json
import os
from pathlib import Path

Path(os.sys.argv[1]).write_text(json.dumps({
    "campaign_id": os.environ["VCFCACHE_CAMPAIGN_ID"],
    "phase": os.environ["VCFCACHE_PHASE"],
    "task_id": int(os.environ["SLURM_ARRAY_TASK_ID"]),
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", os.environ["SLURM_JOB_ID"]),
    "slurm_node": os.environ.get("SLURMD_NODENAME", "unknown"),
}, indent=2, sort_keys=True) + "\n")
PY
mv "$partial_result" "$result_dir"
python3 - "$attempt_dir/success.json" "$result_dir" <<'PY'
import json
import os
import sys
from pathlib import Path

Path(sys.argv[1]).write_text(json.dumps({
    "campaign_id": os.environ["VCFCACHE_CAMPAIGN_ID"],
    "phase": os.environ["VCFCACHE_PHASE"],
    "task_id": int(os.environ["SLURM_ARRAY_TASK_ID"]),
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "result_dir": sys.argv[2],
    "exit_code": 0,
}, indent=2, sort_keys=True) + "\n")
PY
completed=1
