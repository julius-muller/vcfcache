#!/bin/bash
#SBATCH --job-name=vcfcache-adaptive
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
set -euo pipefail

# Arguments are explicit so an immutable deployment can be reused from a fresh
# checkout. Slurm allocation CPUs must equal the requested benchmark CPU count.
deployment=$1
cpus=$2
export PYTHONPATH="$deployment"
export VCFCACHE_ENGINE="$deployment/native/target/release/vcfcache-engine"
cache=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything
output=/results/adaptive-routing/controls-${SLURM_JOB_ID}-cpu${cpus}
args=(--cache "$cache" --output "$output" --cpus "$cpus" --generate-controls --slurm-steps)
case "$cpus" in
  8) args+=(--repeats 3 --warmup --arms legacy original-fused auto regions contigs
            --control-names size-250 span-1 size-75000
            --original-root /results/deployment/vcfcache-fused-c878d86c) ;;
  1) args+=(--repeats 1 --arms legacy auto) ;;
  16) args+=(--repeats 1 --arms legacy auto regions contigs) ;;
  *) echo 'Expected a 1-, 8- or 16-CPU allocation' >&2; exit 2 ;;
esac
"$deployment/.venv/bin/python" -m benchmarks.adaptive_routing "${args[@]}"
