#!/bin/bash
#SBATCH --job-name=native-bcf-clinvar
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=01:00:00
set -euo pipefail
deployment=$1
export PYTHONPATH="$deployment"
root=/mnt/data/tmp/native-veprs-1277
name=LP6005519-DNA_A06.wes
"$deployment/.venv/bin/python" -m benchmarks.native_bcf.diagnose_clinical \
  --cache "$root/full/cache/veprs_clinical" --input "$root/inputs/$name.bcf" \
  --cached "$root/full-runs/$name-native-1/output.bcf" \
  --direct "$root/full-runs/$name-direct-1/output.bcf" \
  --comparison "/results/native-bcf/veprs-clinical-1277/full/$name-direct-1/comparison.json" \
  --output "/results/native-bcf/clinical-diagnostic-${SLURM_JOB_ID}" \
  --work "/mnt/data/tmp/native-clinvar-diagnostic-${SLURM_JOB_ID}"
