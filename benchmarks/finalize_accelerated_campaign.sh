#!/usr/bin/env bash
""":"
Finalize the live accelerated benchmark after external staging.

This controller waits for the original staging controller and assay array,
deploys a fast-forward bundle to the head and workers, retries only missing
assay tasks, and releases the held external jobs. It is intentionally driven by
explicit environment variables so job IDs and commits remain auditable.
":"""

set -euo pipefail

campaign_id=${VCFCACHE_EXTERNAL_CAMPAIGN_ID:?required}
target_commit=${VCFCACHE_TARGET_COMMIT:?required}
bundle=${VCFCACHE_DEPLOYMENT_BUNDLE:?required}
assay_job_id=${VCFCACHE_ASSAY_JOB_ID:?required}
assay_campaign=${VCFCACHE_ASSAY_CAMPAIGN:?required}
external_smoke_job=${VCFCACHE_EXTERNAL_SMOKE_JOB:?required}
external_singlepass_job=${VCFCACHE_EXTERNAL_SINGLEPASS_JOB:?required}
repo_root=${VCFCACHE_REPO_ROOT:-/mnt/data/vcfcache}
results_root=${VCFCACHE_RESULTS_ROOT:-/mnt/data/slurm-results}
identity_file=${VCFCACHE_IDENTITY_FILE:-/home/ubuntu/.ssh/id_ed25519}
staged_marker=$results_root/campaigns/$campaign_id/STAGED
assay_root=$results_root/campaigns/$assay_campaign
worker_ips=(10.0.0.64 10.0.0.75 10.0.0.168 10.0.0.249 10.0.0.159 10.0.0.167)
ssh_options=(-i "$identity_file" -o BatchMode=yes -o StrictHostKeyChecking=accept-new)

while [[ ! -s "$staged_marker" ]]; do
  sleep 30
done

while sudo -u appuser squeue -h -j "$assay_job_id" | grep -q .; do
  sleep 30
done

deploy_local() {
  sudo -u appuser git -C "$repo_root" diff --quiet
  sudo -u appuser git -C "$repo_root" fetch "$bundle" refs/heads/main
  sudo -u appuser git -C "$repo_root" merge --ff-only FETCH_HEAD
  test "$(sudo -u appuser git -C "$repo_root" rev-parse HEAD)" = "$target_commit"
}

deploy_worker() {
  local worker=$1
  ssh "${ssh_options[@]}" "ubuntu@$worker" \
    "sudo -u appuser git -C '$repo_root' diff --quiet && \
     sudo -u appuser git -C '$repo_root' fetch '/results/${bundle#"$results_root"/}' refs/heads/main && \
     sudo -u appuser git -C '$repo_root' merge --ff-only FETCH_HEAD && \
     test \"\$(sudo -u appuser git -C '$repo_root' rev-parse HEAD)\" = '$target_commit'"
}

deploy_local
worker_pids=()
for worker in "${worker_ips[@]}"; do
  deploy_worker "$worker" &
  worker_pids+=("$!")
done
for worker_pid in "${worker_pids[@]}"; do
  wait "$worker_pid"
done

missing=$(python3 - "$assay_root" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
with (root / "manifests/measured.tsv").open(newline="") as handle:
    expected = {int(row["task_id"]) for row in csv.DictReader(handle, delimiter="\t")}
completed = {
    int(path.name.removeprefix("task-"))
    for path in (root / "measured/tasks").glob("task-*")
}
print(",".join(str(value) for value in sorted(expected - completed)))
PY
)

retry_job=""
if [[ -n "$missing" ]]; then
  retry_job=$(sudo -u appuser sbatch --parsable \
    --job-name=assay-semantic-retry \
    --array="$missing%4" \
    --chdir="$repo_root" \
    --output="/results/campaigns/$assay_campaign/logs/retry-%A_%a.out" \
    --export="ALL,VCFCACHE_CAMPAIGN_ID=$assay_campaign,VCFCACHE_PHASE=measured,VCFCACHE_TASK_MANIFEST=/results/campaigns/$assay_campaign/manifests/measured.tsv,VCFCACHE_RESULTS_ROOT=/results,VCFCACHE_REPO_ROOT=$repo_root" \
    "$repo_root/benchmarks/slurm_pair.sh")
  retry_job=${retry_job%%;*}
fi

sudo -u appuser scontrol release "$external_smoke_job" "$external_singlepass_job"

sudo -u appuser python3 - "$results_root/campaigns/$campaign_id/ACCELERATED" \
  "$target_commit" "$retry_job" "$missing" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

Path(sys.argv[1]).write_text(json.dumps({
    "completed_at": datetime.now(timezone.utc).isoformat(),
    "deployed_commit": sys.argv[2],
    "assay_retry_job": sys.argv[3] or None,
    "assay_retry_task_ids": sys.argv[4],
    "external_jobs_released": True,
}, indent=2, sort_keys=True) + "\n")
PY
