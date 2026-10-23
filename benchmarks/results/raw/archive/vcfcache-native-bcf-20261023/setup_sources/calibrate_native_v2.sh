#!/bin/bash
#SBATCH --job-name=native-bcf-calibration
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
export PYTHONPATH="$deployment"
export VCFCACHE_ENGINE="$deployment/native/target/release/vcfcache-engine"
"$deployment/.venv/bin/python" /results/veprs/scripts/calibrate_native_v2.py \
  --cache "$cache" --input "$source" \
  --output "/results/native-bcf/calibration-${SLURM_JOB_ID}" \
  --work-dir "/mnt/data/tmp/native-calibration-${SLURM_JOB_ID}" --cpus 16
