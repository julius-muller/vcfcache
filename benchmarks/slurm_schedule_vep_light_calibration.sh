#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000M
#SBATCH --time=00:20:00

set -euo pipefail

repo_root=${VCFCACHE_REPO_ROOT:?VCFCACHE_REPO_ROOT is required}
campaign_id=${VCFCACHE_VEP_CALIBRATION_CAMPAIGN_ID:?VCFCACHE_VEP_CALIBRATION_CAMPAIGN_ID is required}

export PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}"
python="$repo_root/.venv/bin/python"
runner="$repo_root/benchmarks/run_external_cohort.py"

"$python" "$runner" prepare \
  --campaign-id "$campaign_id" \
  --controller-results /results \
  --worker-results /results \
  --tool vep \
  --calibration-per-cohort 2

"$python" "$runner" submit-chain \
  --campaign-id "$campaign_id" \
  --controller-results /results \
  --worker-results /results \
  --concurrency 3 \
  --smoke-first
