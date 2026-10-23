#!/bin/bash
#SBATCH --job-name=engine-overhead-profile
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28000M
#SBATCH --exclusive
#SBATCH --time=06:00:00
set -euo pipefail
cd /tmp
export PYTHONPATH=$5
/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python \
  -m benchmarks.engine_profile.run group "$1" --group "$2" --attempt "$3" --cells "$4"
