#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --time=00:10:00

set -euo pipefail

repo_root=${VCFCACHE_REPO_ROOT:?VCFCACHE_REPO_ROOT is required}
campaign_root=${VCFCACHE_CAMPAIGN_ROOT:?VCFCACHE_CAMPAIGN_ROOT is required}
phase=${VCFCACHE_PHASE:-measured}

export PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}"
python3 "$repo_root/benchmarks/validate_external_completion.py" \
  --campaign-root "$campaign_root" \
  --phase "$phase" \
  --visibility-attempts 60 \
  --poll-seconds 5
