#!/usr/bin/env bash
# Prepare, stage, and submit the controlled runtime campaign after external WGS.

set -euo pipefail

campaign_id=${VCFCACHE_CAMPAIGN_ID:?required}
target_commit=${VCFCACHE_TARGET_COMMIT:?required}
bundle=${VCFCACHE_DEPLOYMENT_BUNDLE:?required}
prep_bundle=${VCFCACHE_PREP_BUNDLE:?required}
wait_job_id=${VCFCACHE_WAIT_JOB_ID:?required}
prep_host=${VCFCACHE_PREP_HOST:-appuser@10.133.255.21}
prep_root=${VCFCACHE_PREP_ROOT:-/mnt/data/vcfcache_benchmarks}
prep_repo=${VCFCACHE_PREP_REPO:-/mnt/data/home/appuser-projects/vcfcache}
repo_root=${VCFCACHE_REPO_ROOT:-/mnt/data/vcfcache}
data_root=${VCFCACHE_DATA_ROOT:-/mnt/data/vcfcache_benchmarks}
results_root=${VCFCACHE_RESULTS_ROOT:-/mnt/data/slurm-results}
identity_file=${VCFCACHE_IDENTITY_FILE:-/home/ubuntu/.ssh/id_ed25519}
known_hosts=${VCFCACHE_KNOWN_HOSTS:-/home/ubuntu/.ssh/known_hosts}
worker_ips=(10.0.0.64 10.0.0.75 10.0.0.168 10.0.0.249 10.0.0.159 10.0.0.167)
ssh_options=(-i "$identity_file" -o BatchMode=yes -o StrictHostKeyChecking=accept-new -o UserKnownHostsFile="$known_hosts")
ssh_command="ssh -i $identity_file -o BatchMode=yes -o StrictHostKeyChecking=accept-new -o UserKnownHostsFile=$known_hosts"
controlled_root=$data_root/controlled_runtime
campaign_root=$results_root/campaigns/$campaign_id

# The external-preparation command writes READY only after normalization, QC,
# relatedness screening, and all three cohort caches are complete. Updating its
# checkout before that atomic marker would make the live run non-reproducible.
while ! ssh "${ssh_options[@]}" "$prep_host" \
  "test -s '$prep_root/external_wgs/READY.json'"; do
  sleep 60
done

ssh "${ssh_options[@]}" "$prep_host" \
  "git -C '$prep_repo' diff --quiet && \
   git -C '$prep_repo' fetch '$prep_bundle' refs/heads/main && \
   git -C '$prep_repo' merge --ff-only FETCH_HEAD && \
   test \"\$(git -C '$prep_repo' rev-parse HEAD)\" = '$target_commit' && \
   cd '$prep_repo' && .venv/bin/python benchmarks/prepare_controlled_runtime.py \
     >'$prep_root/controlled_runtime.prepare.log' 2>&1"
ssh "${ssh_options[@]}" "$prep_host" \
  "test -s '$prep_root/controlled_runtime/READY.json'"

# Never change the worker checkout underneath an external benchmark process.
while sudo -u appuser squeue -h -j "$wait_job_id" | grep -q .; do
  sleep 60
done

mkdir -p "$controlled_root"
rsync -aH --delete --partial -e "$ssh_command" \
  "$prep_host:$prep_root/controlled_runtime/" "$controlled_root/"
chown -R appuser:appgroup "$controlled_root"

stage_worker() {
  local worker=$1
  rsync -aH --delete --partial --rsync-path="sudo rsync" -e "$ssh_command" \
    "$controlled_root/" "ubuntu@$worker:$controlled_root/"
  ssh "${ssh_options[@]}" "ubuntu@$worker" \
    "sudo chown -R appuser:appgroup '$controlled_root'"
}

worker_pids=()
for worker in "${worker_ips[@]}"; do
  stage_worker "$worker" &
  worker_pids+=("$!")
done
for worker_pid in "${worker_pids[@]}"; do
  wait "$worker_pid"
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

sudo -u appuser "$repo_root/.venv/bin/python" \
  "$repo_root/benchmarks/run_controlled_cohort.py" prepare \
  --campaign-id "$campaign_id"
submission=$(sudo -u appuser "$repo_root/.venv/bin/python" \
  "$repo_root/benchmarks/run_controlled_cohort.py" submit \
  --campaign-id "$campaign_id" --concurrency 6)

sudo -u appuser python3 - "$campaign_root/CONTROLLED_SUBMITTED" \
  "$target_commit" "$wait_job_id" "$submission" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

Path(sys.argv[1]).write_text(json.dumps({
    "submitted_at": datetime.now(timezone.utc).isoformat(),
    "commit": sys.argv[2],
    "waited_for_external_job": sys.argv[3],
    "submission": json.loads(sys.argv[4]),
}, indent=2, sort_keys=True) + "\n")
PY
