#!/bin/bash
#SBATCH --job-name=native-bcf-veprs-cache
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=04:00:00
set -euo pipefail
deployment=$1
export PYTHONPATH="$deployment"
"$deployment/.venv/bin/python" "$2" --output "$3" --transcripts /results/veprs/json_cache/GRCh38_native_20261023
