#!/bin/bash
#SBATCH --job-name=branch-eval
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
# Time limit is supplied at submission; annotators differ by two orders of magnitude.
#
# One annotator, one sample, nested input scales, every engine generation.
# `direct` runs first so that uncached annotation, not another engine, is the
# semantic reference every cached arm is audited against.
set -euo pipefail
annotator=$1       # vep | fastvep | veprs
scales=$2          # space-separated subset of: single minimal panel wes wgs
arms=$3            # space-separated subset of the driver's arm names
cpus=${4:-16}

head_root=/results/deployment/vcfcache-eval-head-63757a07
main_root=/results/deployment/vcfcache-eval-main-136bafa9
fused_root=/results/deployment/vcfcache-fused-c878d86c
runner=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python

case "$annotator" in
  vep)
    origin=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh38_joint_af001
    name=vep115.2_everything ;;
  fastvep)
    origin=/results/richness/caches/fastvep_clinical_grch38_af001
    name=fastvep_clinical_grch38_af001 ;;
  veprs)
    origin=/results/native-bcf/veprs-caches-20261023b/veprs_core
    name=veprs_core ;;
  veprs-clinical)
    # The heavy vep-rs rung: genome-wide gnomAD plus ClinVar, ~40 us/variant
    # against ~12 us for core. Its cached output is known to disagree with fresh
    # annotation at 22/82,714 records in ClinVar-derived CSQ values; the audit
    # records that rather than the run being abandoned over it.
    origin=/results/native-bcf/veprs-caches-20261023b/veprs_clinical
    name=veprs_clinical ;;
  *) echo "unknown annotator: $annotator" >&2; exit 2 ;;
esac

work=/mnt/data/tmp/branch-eval-${annotator}-${SLURM_JOB_ID}
output=/results/branch-eval/${annotator}-${SLURM_JOB_ID}
mkdir -p "$work/cache/$name" "$work/workflow" "$work/blueprint" "$output"

cp "$origin/workflow/init.yaml" "$work/workflow/"
cp "$origin/cache/$name/"{vcfcache_annotated.bcf,vcfcache_annotated.bcf.csi,annotation.yaml,params.snapshot.yaml} \
   "$work/cache/$name/"
for extra in build_provenance.json; do
  if [ -f "$origin/cache/$name/$extra" ]; then cp "$origin/cache/$name/$extra" "$work/cache/$name/"; fi
done
# The vep-rs pilot recipe carries a params-only key that the annotation schema
# rejects; repair the staged copy only, never the source cache.
if grep -q '^tool_version_command:' "$work/cache/$name/annotation.yaml"; then
  PYTHONPATH="$head_root" "$runner" /results/veprs/scripts/repair_native_recipe.py "$work/cache/$name"
fi

inputs=()
for scale in $scales; do inputs+=("/results/branch-eval/inputs/$scale.bcf"); done

sha256sum "$work/cache/$name/"{vcfcache_annotated.bcf,annotation.yaml,params.snapshot.yaml} \
  "$head_root/native/target/release/vcfcache-engine" > "$output/SOURCE_SHA256SUMS"
{ echo "annotator=$annotator"; echo "scales=$scales"; echo "arms=$arms"; echo "cpus=$cpus";
  echo "node=$(hostname)"; echo "job=${SLURM_JOB_ID}"; } > "$output/CONDITIONS.txt"

# `python -c` puts the working directory ahead of PYTHONPATH on sys.path, so a
# working directory that contains a vcfcache package would override the
# per-arm deployment the driver selects. Run from a directory that has none,
# and prove each deployment still resolves to itself before timing anything.
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
  --cache "$work/cache/$name" \
  --input "${inputs[@]}" \
  --output "$output/benchmark" \
  --work-dir "$work/runs" \
  --cpus "$cpus" --repeats 1 --slurm-steps \
  --arms $arms \
  --original-root "$fused_root" \
  --legacy-root "$main_root"
