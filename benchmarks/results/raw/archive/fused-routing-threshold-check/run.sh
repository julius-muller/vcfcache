#!/bin/bash
# Two record counts, same condition, so per-variant cost is a fitted slope for
# each engine rather than one wall time. Work dir is on NVMe, not tmpfs: only
# ~6 GB of RAM is free and the fused engine deliberately trades temp disk for
# fewer full-file scans, so a disk-backed run is the fairer test of that trade.
set -euo pipefail
cd /home/j380r/projects/vcfcache
export VCFCACHE_ENGINE="$PWD/native/target/release/vcfcache-engine"
BASE=/home/j380r/.cache/vcfcache-fused-threshold-check
mkdir -p "$BASE"
for N in 200000 800000; do
  echo "### records=$N"
  .venv/bin/python -m benchmarks.fused_engine \
    --output "$BASE/n$N" --work-dir "$BASE/work-n$N" \
    --records "$N" --hit-rates 0.9 --payloads 1024 --cpus 8
done
echo "### summaries"
for N in 200000 800000; do
  echo "--- n$N"; cat "$BASE/n$N/summary.tsv" 2>/dev/null || echo "(none)"
done
