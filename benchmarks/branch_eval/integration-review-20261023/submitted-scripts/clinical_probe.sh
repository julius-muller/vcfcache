#!/bin/bash
#SBATCH --job-name=integration-clinical-probe
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=01:00:00
set -euo pipefail
python=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
"$python" "$1" --output "/results/branch-integration/clinical-${SLURM_JOB_ID}" \
  --work "/mnt/data/tmp/integration-clinical-${SLURM_JOB_ID}"
