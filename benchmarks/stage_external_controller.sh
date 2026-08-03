#!/usr/bin/env bash
# Stage a prepared external campaign from the preparation VM onto the cluster.

set -euo pipefail

campaign_id=${VCFCACHE_CAMPAIGN_ID:?VCFCACHE_CAMPAIGN_ID is required}
# A Slurm gate can preserve an existing dependency chain, but it must not be
# required: waiting for preparation does not need to reserve an entire worker.
gate_job_id=${VCFCACHE_GATE_JOB_ID:-}
smoke_job_id=${VCFCACHE_SMOKE_JOB_ID:?VCFCACHE_SMOKE_JOB_ID is required}
warmup_job_id=${VCFCACHE_WARMUP_JOB_ID:?VCFCACHE_WARMUP_JOB_ID is required}
measured_job_id=${VCFCACHE_MEASURED_JOB_ID:?VCFCACHE_MEASURED_JOB_ID is required}
target_commit=${VCFCACHE_TARGET_COMMIT:?VCFCACHE_TARGET_COMMIT is required}
base_commit=${VCFCACHE_BASE_COMMIT:-32a760323fc5}
prep_host=${VCFCACHE_PREP_HOST:-appuser@10.133.255.21}
prep_root=${VCFCACHE_PREP_ROOT:-/mnt/data/vcfcache_benchmarks}
repo_root=${VCFCACHE_REPO_ROOT:-/mnt/data/vcfcache}
data_root=${VCFCACHE_DATA_ROOT:-/mnt/data/vcfcache_benchmarks}
vep_root=${VCFCACHE_VEP_ROOT:-/mnt/data/apps/ensembl-vep/115}
results_root=${VCFCACHE_RESULTS_ROOT:-/mnt/data/slurm-results}
identity_file=${VCFCACHE_IDENTITY_FILE:-/home/ubuntu/.ssh/id_ed25519}
known_hosts=${VCFCACHE_KNOWN_HOSTS:-/home/ubuntu/.ssh/known_hosts}
worker_ips=(10.0.0.64 10.0.0.75 10.0.0.168 10.0.0.249 10.0.0.159 10.0.0.167)
ssh_command="ssh -i $identity_file -o BatchMode=yes -o StrictHostKeyChecking=accept-new -o UserKnownHostsFile=$known_hosts"
external_root=$data_root/external_wgs
campaign_root=$results_root/campaigns/$campaign_id
deferred_root=$results_root/deferred/$campaign_id

if [[ -n "$gate_job_id" ]]; then
  while :; do
    state=$(sudo -u appuser squeue -h -j "$gate_job_id" -o %T | head -1)
    if [[ "$state" == RUNNING ]]; then
      break
    fi
    if [[ -z "$state" ]]; then
      echo "Staging gate $gate_job_id left the queue before running" >&2
      exit 1
    fi
    sleep 30
  done
fi

while ! $ssh_command "$prep_host" "test -s $prep_root/external_wgs/READY.json"; do
  sleep 60
done

mkdir -p "$external_root/deployment"
rsync -a --partial -e "$ssh_command" \
  "$prep_host:$prep_root/external_wgs/deployment/vcfcache-external.bundle" \
  "$external_root/deployment/vcfcache-external.bundle"
for relative in samples qc manifests cohort_caches; do
  mkdir -p "$external_root/$relative"
  rsync -a --partial -e "$ssh_command" \
    "$prep_host:$prep_root/external_wgs/$relative/" \
    "$external_root/$relative/"
done
rsync -a -e "$ssh_command" \
  "$prep_host:$prep_root/external_wgs/READY.json" "$external_root/READY.json"
for cache in gnomad_v4.1_GRCh37_joint_af010 gnomad_v4.1_GRCh37_joint_af001; do
  mkdir -p "$data_root/bundled_zenodo_caches/$cache"
  rsync -a --partial -e "$ssh_command" \
    "$prep_host:$prep_root/bundled_zenodo_caches/$cache/" \
    "$data_root/bundled_zenodo_caches/$cache/"
done
chown -R appuser:appgroup "$external_root" "$data_root/bundled_zenodo_caches"
test -s "$vep_root/cachedir/homo_sapiens/115_GRCh37/info.txt"
test -s "$vep_root/cachedir/homo_sapiens/115_GRCh37/1/all_vars.gz"

current=$(sudo -u appuser git -C "$repo_root" rev-parse --short=12 HEAD)
if [[ "$current" != "$target_commit" ]]; then
  [[ "$current" == "$base_commit" ]]
  sudo -u appuser git -C "$repo_root" fetch \
    "$external_root/deployment/vcfcache-external.bundle" refs/heads/main
  sudo -u appuser git -C "$repo_root" merge --ff-only FETCH_HEAD
fi

stage_worker() {
  local worker=$1
  local remote="ubuntu@$worker"
  for relative in samples qc manifests cohort_caches; do
    rsync -a --partial --rsync-path="sudo rsync" -e "$ssh_command" \
      "$external_root/$relative/" "$remote:$external_root/$relative/"
  done
  rsync -a --rsync-path="sudo rsync" -e "$ssh_command" \
    "$external_root/READY.json" "$remote:$external_root/READY.json"
  rsync -a --rsync-path="sudo rsync" -e "$ssh_command" \
    "$external_root/deployment/vcfcache-external.bundle" \
    "$remote:$external_root/deployment/vcfcache-external.bundle"
  for cache in gnomad_v4.1_GRCh37_joint_af010 gnomad_v4.1_GRCh37_joint_af001; do
    rsync -a --partial --rsync-path="sudo rsync" -e "$ssh_command" \
      "$data_root/bundled_zenodo_caches/$cache/" \
      "$remote:$data_root/bundled_zenodo_caches/$cache/"
  done
  rsync -a --partial --rsync-path="sudo rsync" -e "$ssh_command" \
    "$vep_root/cachedir/homo_sapiens/115_GRCh37/" \
    "$remote:$vep_root/cachedir/homo_sapiens/115_GRCh37/"
  $ssh_command "$remote" "sudo chown -R appuser:appgroup '$external_root' '$data_root/bundled_zenodo_caches' '$vep_root/cachedir/homo_sapiens/115_GRCh37'; test -s '$vep_root/cachedir/homo_sapiens/115_GRCh37/info.txt'; test -s '$vep_root/cachedir/homo_sapiens/115_GRCh37/1/all_vars.gz'; current=\$(sudo -u appuser git -C '$repo_root' rev-parse --short=12 HEAD); if [[ \"\$current\" != '$target_commit' ]]; then [[ \"\$current\" == '$base_commit' ]]; sudo -u appuser git -C '$repo_root' fetch '$external_root/deployment/vcfcache-external.bundle' refs/heads/main; sudo -u appuser git -C '$repo_root' merge --ff-only FETCH_HEAD; fi; test \"\$(sudo -u appuser git -C '$repo_root' rev-parse --short=12 HEAD)\" = '$target_commit'"
}

worker_pids=()
for worker in "${worker_ips[@]}"; do
  stage_worker "$worker" &
  worker_pids+=("$!")
done
for worker_pid in "${worker_pids[@]}"; do
  wait "$worker_pid"
done

sudo -u appuser "$repo_root/.venv/bin/python" \
  "$repo_root/benchmarks/run_external_cohort.py" prepare \
  --campaign-id "$campaign_id"
cp "$deferred_root/submission.json" "$campaign_root/submission.json"
chown appuser:appgroup "$campaign_root/submission.json"
python3 - "$campaign_root/STAGED" <<PY
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

Path(sys.argv[1]).write_text(json.dumps({
    "campaign_id": "$campaign_id",
    "created_at": datetime.now(timezone.utc).isoformat(),
    "commit": "$target_commit",
    "gate_job_id": "$gate_job_id",
    "smoke_job_id": "$smoke_job_id",
    "warmup_job_id": "$warmup_job_id",
    "measured_job_id": "$measured_job_id",
    "workers": ${#worker_ips[@]},
}, indent=2, sort_keys=True) + "\n")
PY
chown appuser:appgroup "$campaign_root/STAGED"
