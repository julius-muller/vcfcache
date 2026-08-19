#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000M
#SBATCH --time=00:20:00

set -euo pipefail

repo_root=${VCFCACHE_REPO_ROOT:?VCFCACHE_REPO_ROOT is required}
campaign_id=${VCFCACHE_FASTVEP_CAMPAIGN_ID:?VCFCACHE_FASTVEP_CAMPAIGN_ID is required}
strategy_manifest=/mnt/data/vcfcache_benchmarks/fastvep_external_caches/fastvep_strategies.json

export PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}"
/mnt/data/vcfcache/.venv/bin/python "$repo_root/benchmarks/run_external_cohort.py" \
  prepare \
  --campaign-id "$campaign_id" \
  --controller-results /results \
  --worker-results /results \
  --tool fastvep \
  --strategy-manifest "$strategy_manifest"

/mnt/data/vcfcache/.venv/bin/python "$repo_root/benchmarks/run_external_cohort.py" \
  submit-chain \
  --campaign-id "$campaign_id" \
  --controller-results /results \
  --worker-results /results \
  --concurrency 6 \
  --smoke-first
