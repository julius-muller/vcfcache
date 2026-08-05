#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --exclusive
#SBATCH --time=02:00:00

set -euo pipefail

repo_root=${VCFCACHE_REPO_ROOT:?VCFCACHE_REPO_ROOT is required}
campaign_root=${VCFCACHE_CAMPAIGN_ROOT:?VCFCACHE_CAMPAIGN_ROOT is required}
phase=${VCFCACHE_PHASE:-warmup}
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}

export LC_ALL=C
export LANG=C

srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=8 \
  "$repo_root/.venv/bin/python" \
  "$repo_root/benchmarks/recover_external_attempts.py" \
  --campaign-root "$campaign_root" \
  --phase "$phase" \
  --task-id "$task_id"
