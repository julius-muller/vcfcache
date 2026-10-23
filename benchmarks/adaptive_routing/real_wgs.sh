#!/bin/bash
#SBATCH --job-name=vcfcache-adaptive-wgs
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=03:00:00
set -euo pipefail
deployment=$1
sample=$2
export PYTHONPATH="$deployment"
export VCFCACHE_ENGINE="$deployment/native/target/release/vcfcache-engine"
original=/results/richness/caches/fastvep_clinical_grch38_af001
name=fastvep_clinical_grch38_af001
work=/mnt/data/tmp/adaptive-real-${SLURM_JOB_ID}
mkdir -p "$work/base/blueprint" "$work/base/workflow" "$work/base/cache/$name"
# Stage the unchanged source BCF on local storage for BOTH candidate engines.
cp "$original/workflow/init.yaml" "$work/base/workflow/"
cp "$original/cache/$name/"{vcfcache_annotated.bcf,vcfcache_annotated.bcf.csi,annotation.yaml,params.snapshot.yaml} "$work/base/cache/$name/"
cache=$work/base/cache/$name
source=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/sgdp/${sample}.GRCh38.small_variants.vcf.gz
test -s "$source"
output=/results/adaptive-routing/wgs-${SLURM_JOB_ID}-${sample}
mkdir -p "$output"
sha256sum "$cache/"{vcfcache_annotated.bcf,annotation.yaml,params.snapshot.yaml} > "$output/SHA256SUMS"
"$VCFCACHE_ENGINE" store-build --cache "$cache/vcfcache_annotated.bcf" --output "$work/store"
cp "$work/store/manifest.json" "$output/store_manifest.json"
"$deployment/.venv/bin/python" -m benchmarks.adaptive_routing --cache "$cache" --input "$source" \
  --output "$output/benchmark" --cpus 16 --repeats 1 --arms legacy auto store \
  --store "$work/store" --work-dir "$work/runs" --slurm-steps
# Preserve local store/work for diagnosis. Cleanup is explicit after collection.
