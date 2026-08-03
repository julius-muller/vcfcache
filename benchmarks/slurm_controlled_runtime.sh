#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --exclusive
#SBATCH --time=12:00:00

set -euo pipefail

repo_root=${VCFCACHE_REPO_ROOT:?required}
task_manifest=${VCFCACHE_TASK_MANIFEST:?required}
results_root=${VCFCACHE_RESULTS_ROOT:-/results}
campaign_id=${VCFCACHE_CAMPAIGN_ID:?required}
phase=${VCFCACHE_PHASE:?required}
task_id=${SLURM_ARRAY_TASK_ID:?required}
run_root=/mnt/data/tmp/vcfcache_controlled_runs/$campaign_id/$phase/task-$task_id
phase_root=$results_root/campaigns/$campaign_id/$phase
result_dir=$phase_root/tasks/task-$task_id
attempt_dir=$phase_root/attempts/task-$task_id/job-$SLURM_JOB_ID
partial_result=$result_dir.partial-$SLURM_JOB_ID
completed=0

archive_attempt() {
  local status=$?
  if ((completed)); then return 0; fi
  mkdir -p "$attempt_dir"
  if [[ -d "$run_root" ]]; then rsync -aH "$run_root/" "$attempt_dir/run/"; fi
  python3 - "$attempt_dir/failure.json" "$status" <<'PY'
import json, os, sys
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps({
    "campaign_id": os.environ["VCFCACHE_CAMPAIGN_ID"],
    "phase": os.environ["VCFCACHE_PHASE"],
    "task_id": int(os.environ["SLURM_ARRAY_TASK_ID"]),
    "job_id": os.environ["SLURM_JOB_ID"],
    "node": os.environ.get("SLURMD_NODENAME", "unknown"),
    "exit_code": int(sys.argv[2]),
}, indent=2, sort_keys=True) + "\n")
PY
  return "$status"
}
trap archive_attempt EXIT

if [[ -d "$result_dir" ]]; then completed=1; exit 0; fi
rm -rf -- "$run_root"
mkdir -p "$run_root" "$attempt_dir"
export LC_ALL=C LANG=C VCFCACHE_REQUIRE_CGROUP_METRICS=1
srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=8 \
  "$repo_root/.venv/bin/python" "$repo_root/benchmarks/run_controlled_task.py" \
  --task-manifest "$task_manifest" --task-id "$task_id" --run-root "$run_root"
test -s "$run_root/controlled_summary.json"
test ! -e "$partial_result"
mkdir "$partial_result"
rsync -aH "$run_root/" "$partial_result/"
mv "$partial_result" "$result_dir"
python3 - "$attempt_dir/success.json" "$result_dir" <<'PY'
import json, os, sys
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps({
    "task_id": int(os.environ["SLURM_ARRAY_TASK_ID"]),
    "job_id": os.environ["SLURM_JOB_ID"],
    "result_dir": sys.argv[2],
    "exit_code": 0,
}, indent=2, sort_keys=True) + "\n")
PY
completed=1
