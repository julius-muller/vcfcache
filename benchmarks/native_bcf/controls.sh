#!/bin/bash
#SBATCH --job-name=native-bcf-controls
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=01:00:00
set -euo pipefail
deployment=$1
export PYTHONPATH="$deployment"
export VCFCACHE_ENGINE="$deployment/native/target/release/vcfcache-engine"
name=fastvep_clinical_grch38_af001
original=/results/richness/caches/$name
work=/mnt/data/tmp/native-controls-${SLURM_JOB_ID}
output=/results/native-bcf/controls-${SLURM_JOB_ID}
mkdir -p "$work/base/blueprint" "$work/base/workflow" "$work/base/cache/$name" "$work/inputs"
cp "$original/workflow/init.yaml" "$work/base/workflow/"
cp "$original/cache/$name/"{vcfcache_annotated.bcf,vcfcache_annotated.bcf.csi,annotation.yaml,params.snapshot.yaml} "$work/base/cache/$name/"
bed=/results/richness/regions/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed
source=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/sgdp/LP6005519-DNA_A06.GRCh38.small_variants.vcf.gz
bcftools view -R "$bed" "$source" -Ob -W -o "$work/wes.bcf"
bcftools norm -m- "$work/wes.bcf" -Ob -W -o "$work/normalized.bcf"
bcftools annotate -a "$work/base/cache/$name/vcfcache_annotated.bcf" "$work/normalized.bcf" -c INFO -Ou |
  bcftools view -i 'INFO/CSQ!=""' -Ou |
  bcftools annotate -x INFO -Ob -W -o "$work/inputs/all-hit-control.bcf"
"$deployment/.venv/bin/python" -m benchmarks.adaptive_routing \
  --cache "$work/base/cache/$name" --input "$work/inputs/all-hit-control.bcf" \
  --output "$output/benchmark" --work-dir "$work/runs" --cpus 16 \
  --repeats 3 --paired-order --arms auto native --slurm-steps
# WGS buffer-budget compatibility smoke test. It is not a new WGS speedup claim.
"$deployment/.venv/bin/python" -m benchmarks.adaptive_routing \
  --cache "$work/base/cache/$name" --input "$source" \
  --output "$output/wgs" --work-dir "$work/wgs-runs" --cpus 16 \
  --repeats 1 --arms auto native --slurm-steps
