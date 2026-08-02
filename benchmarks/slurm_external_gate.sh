#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --time=24:00:00

set -euo pipefail

gate_marker=${VCFCACHE_GATE_MARKER:?VCFCACHE_GATE_MARKER is required}
deadline=$((SECONDS + 23 * 60 * 60))

while [[ ! -s "$gate_marker" ]]; do
  if ((SECONDS >= deadline)); then
    echo "Timed out waiting for external benchmark staging: $gate_marker" >&2
    exit 1
  fi
  sleep 30
done

echo "External benchmark staging gate opened: $gate_marker"
