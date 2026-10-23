#!/bin/bash
#SBATCH --job-name=vcfcache-adaptive-io
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --time=00:30:00
set -euo pipefail
deployment=$1
export PYTHONPATH="$deployment"
export VCFCACHE_ENGINE="$deployment/native/target/release/vcfcache-engine"
cache=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything
input=/results/adaptive-routing/controls-1258-cpu16/inputs/size-250.bcf
params=/results/adaptive-routing/controls-1257-cpu8/params.yaml
output=/results/adaptive-routing/io-${SLURM_JOB_ID}
mkdir "$output"
for engine in legacy fused; do
  mkdir "$output/$engine"
  srun --exclusive --exact --cpu-bind=cores -n1 -c8 "$deployment/.venv/bin/python" \
    "$deployment/benchmarks/adaptive_routing/trace_cache_io.py" \
    --cache-bcf "$cache/vcfcache_annotated.bcf" --trace "$output/$engine.trace" \
    --summary "$output/$engine.io.json" -- \
    "$deployment/.venv/bin/python" -m vcfcache.cli annotate \
    -a "$cache" -i "$input" -o "$output/$engine/output.bcf" -y "$params" \
    --stats-dir "$output/$engine" --statistics light --engine "$engine" --debug
done
# Traced wall times are deliberately not used for performance comparisons.
