#!/bin/bash
#SBATCH --job-name=branch-eval-heavy
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#
# The heaviest VEP pipelines already measured: VEP --everything plus a
# SyntheticDelay plugin, on KPGP-00319 WGS. The uncached arm of this campaign
# costs 22.7 h at the 10 ms rung, so it is reused from the archive rather than
# rerun; only the engine arms run here.
#
# CPUs are 8, not the 16 used by the assay matrix, because the archived campaign
# allocated 8 and set vep_forks=8. Matching it is what makes the archived
# uncached number comparable; the driver's parameter overrides then reproduce
# the campaign's params.yaml exactly.
#
# All arms run in one job so that every engine comparison is within-job. The
# first arm is the semantic reference for the others, which makes this a genuine
# cross-engine equivalence check even though no fresh uncached arm exists.
set -euo pipefail
pipeline=$1        # e.g. everything_delay_10000us
arms=$2            # driver arm names: the adaptive engine is called 'auto'

cpus=${3:-8}

head_root=/results/deployment/vcfcache-eval-head-63757a07
main_root=/results/deployment/vcfcache-eval-main-136bafa9
fused_root=/results/deployment/vcfcache-fused-c878d86c
runner=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
origin=/mnt/data/vcfcache_benchmarks/wgs_pipeline_spectrum_af001
input=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp/KPGP-00319.GRCh38.small_variants.vcf.gz

work=/mnt/data/tmp/branch-eval-heavy-${pipeline}-${SLURM_JOB_ID}
output=/results/branch-eval/heavy-${pipeline}-${SLURM_JOB_ID}
mkdir -p "$work/cache/$pipeline" "$work/workflow" "$work/blueprint" "$output"

cp "$origin/workflow/init.yaml" "$work/workflow/"
# The campaign cache BCF is a symlink into the bundled Zenodo cache; copy
# dereferences it so the engine reads a local file, as the campaign did.
cp -L "$origin/cache/$pipeline/"{vcfcache_annotated.bcf,vcfcache_annotated.bcf.csi,annotation.yaml,params.snapshot.yaml} \
   "$work/cache/$pipeline/"

sha256sum "$work/cache/$pipeline/"{annotation.yaml,params.snapshot.yaml} \
  "$head_root/native/target/release/vcfcache-engine" > "$output/SOURCE_SHA256SUMS"
{ echo "pipeline=$pipeline"; echo "arms=$arms"; echo "cpus=$cpus"; echo "input=$input";
  echo "node=$(hostname)"; echo "job=${SLURM_JOB_ID}";
  echo "archived_uncached_campaign=wgs-pipeline-spectrum-af001-cd87186-kpgp00319"; } > "$output/CONDITIONS.txt"

cd "$work"
for pair in "head:$head_root" "main:$main_root" "fused:$fused_root"; do
  label=${pair%%:*}
  deployment=${pair#*:}
  resolved=$(PYTHONPATH="$deployment" "$runner" -c 'import vcfcache; print(vcfcache.__file__)')
  case "$resolved" in
    "$deployment"/*) printf 'arm source %s -> %s\n' "$label" "$resolved" ;;
    *) echo "deployment $label resolved to $resolved, not $deployment" >&2; exit 3 ;;
  esac
done > "$output/ARM_SOURCES.txt"
cat "$output/ARM_SOURCES.txt"

PYTHONPATH="$head_root" "$runner" -m benchmarks.adaptive_routing \
  --cache "$work/cache/$pipeline" \
  --input "$input" \
  --output "$output/benchmark" \
  --work-dir "$work/runs" \
  --cpus "$cpus" --repeats 1 --slurm-steps \
  --arms $arms \
  --original-root "$fused_root" \
  --legacy-root "$main_root"
