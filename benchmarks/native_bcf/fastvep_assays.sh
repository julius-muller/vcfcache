#!/bin/bash
#SBATCH --job-name=native-bcf-fastvep
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=02:00:00
set -euo pipefail
deployment=$1
export PYTHONPATH="$deployment"
export VCFCACHE_ENGINE="$deployment/native/target/release/vcfcache-engine"
original=/results/richness/caches/fastvep_clinical_grch38_af001
name=fastvep_clinical_grch38_af001
work=/mnt/data/tmp/native-fastvep-${SLURM_JOB_ID}
output=/results/native-bcf/fastvep-${SLURM_JOB_ID}
mkdir -p "$work/full/blueprint" "$work/full/workflow" "$work/full/cache/$name" "$work/inputs" "$output"
cp "$original/workflow/init.yaml" "$work/full/workflow/"
cp "$original/cache/$name/"{vcfcache_annotated.bcf,vcfcache_annotated.bcf.csi,annotation.yaml,params.snapshot.yaml} "$work/full/cache/$name/"
bed=/results/richness/regions/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed
for sample in LP6005519-DNA_A06 LP6005592-DNA_C05; do
  source=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/sgdp/${sample}.GRCh38.small_variants.vcf.gz
  bcftools view -R "$bed" "$source" -Ob -W -o "$work/inputs/${sample}.wes.bcf"
done
sha256sum "$work/full/cache/$name/"{vcfcache_annotated.bcf,annotation.yaml,params.snapshot.yaml} "$bed" > "$output/SOURCE_SHA256SUMS"
"$deployment/.venv/bin/python" -m benchmarks.adaptive_routing \
  --cache "$work/full/cache/$name" --input "$work/inputs/"*.wes.bcf \
  --output "$output/full" --work-dir "$work/full-runs" --cpus 16 --repeats 1 \
  --arms auto native direct --slurm-steps
# Only construct/measure the projection after recording the generic path.
mkdir -p "$work/projected/blueprint" "$work/projected/workflow" "$work/projected/cache/$name"
cp "$original/workflow/init.yaml" "$work/projected/workflow/"
cp "$original/cache/$name/"{annotation.yaml,params.snapshot.yaml} "$work/projected/cache/$name/"
start=$(date +%s.%N)
bcftools view -R "$bed" --regions-overlap 0 "$work/full/cache/$name/vcfcache_annotated.bcf" \
  -Ob -W --threads 8 -o "$work/projected/cache/$name/vcfcache_annotated.bcf"
end=$(date +%s.%N)
awk -v start="$start" -v end="$end" 'BEGIN {print end-start}' > "$output/projection_build_seconds.txt"
sha256sum "$work/projected/cache/$name/"{vcfcache_annotated.bcf,annotation.yaml} > "$output/PROJECTION_SHA256SUMS"
du -b "$work/projected/cache/$name/"{vcfcache_annotated.bcf,vcfcache_annotated.bcf.csi} > "$output/projection_bytes.tsv"
"$deployment/.venv/bin/python" -m benchmarks.adaptive_routing \
  --cache "$work/projected/cache/$name" --input "$work/inputs/"*.wes.bcf \
  --output "$output/projected" --work-dir "$work/projected-runs" --cpus 16 --repeats 1 \
  --arms auto native direct --slurm-steps
