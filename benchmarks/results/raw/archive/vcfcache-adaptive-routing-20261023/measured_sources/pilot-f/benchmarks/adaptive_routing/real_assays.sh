#!/bin/bash
#SBATCH --job-name=vcfcache-adaptive-assays
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
work=/mnt/data/tmp/adaptive-assays-${SLURM_JOB_ID}
mkdir -p "$work/base/blueprint" "$work/base/workflow" "$work/base/cache/$name" "$work/inputs"
cp "$original/workflow/init.yaml" "$work/base/workflow/"
cp "$original/cache/$name/"{vcfcache_annotated.bcf,vcfcache_annotated.bcf.csi,annotation.yaml,params.snapshot.yaml} "$work/base/cache/$name/"
source=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/sgdp/LP6005519-DNA_A06.GRCh38.small_variants.vcf.gz
bcftools view -R /results/richness/regions/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed \
  "$source" -Ob -W -o "$work/inputs/panel.bcf"
bcftools view -R /results/richness/regions/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed \
  "$source" -Ob -W -o "$work/inputs/wes.bcf"
"$deployment/.venv/bin/python" -m benchmarks.adaptive_routing \
  --cache "$work/base/cache/$name" --input "$work/inputs/panel.bcf" "$work/inputs/wes.bcf" \
  --output "/results/adaptive-routing/assays-${SLURM_JOB_ID}/benchmark" --work-dir "$work/runs" \
  --cpus 16 --repeats 1 --arms legacy auto direct --slurm-steps
