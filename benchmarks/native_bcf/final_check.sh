#!/bin/bash
#SBATCH --job-name=native-bcf-final
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=01:00:00
set -euo pipefail
deployment=$1
cache=$2
source=$3
reference=$4
export PYTHONPATH="$deployment"
export VCFCACHE_ENGINE="$deployment/native/target/release/vcfcache-engine"
"$deployment/.venv/bin/python" -m benchmarks.adaptive_routing \
  --cache "$cache" --input "$source" --reference-root "$reference" \
  --output "/results/native-bcf/final-${SLURM_JOB_ID}" \
  --work-dir "/mnt/data/tmp/native-final-${SLURM_JOB_ID}" --cpus 16 \
  --repeats 1 --arms native --slurm-steps
