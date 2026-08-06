#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --exclusive
#SBATCH --time=12:00:00

set -euo pipefail

repo_root=${VCFCACHE_REPO_ROOT:?VCFCACHE_REPO_ROOT is required}
shared_root=${VCFCACHE_FASTVEP_SHARED:?VCFCACHE_FASTVEP_SHARED is required}
vep_strategies=${VCFCACHE_VEP_STRATEGIES:?VCFCACHE_VEP_STRATEGIES is required}
local_root=/mnt/data/apps/fastvep/0.3.0-publication
cache_root=/mnt/data/vcfcache_benchmarks/fastvep_external_caches

export LC_ALL=C
export LANG=C
export PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p /mnt/data/tmp/fastvep-publication

srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=8 \
  /mnt/data/vcfcache/.venv/bin/python \
  "$repo_root/benchmarks/prepare_fastvep_node.py" \
  --root "$local_root" \
  --cache-root "$cache_root" \
  --fastvep-binary "$shared_root/fastvep" \
  --vep-strategies "$vep_strategies" \
  --vcfcache "$repo_root/.venv/bin/vcfcache" \
  --threads 8
