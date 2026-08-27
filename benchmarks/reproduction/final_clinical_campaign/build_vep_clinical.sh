#!/bin/bash
# Build one immutable matched-clinical VEP cache per assembly.
#SBATCH --job-name=vep-clin-build
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000M
#SBATCH --time=24:00:00
#SBATCH --exclusive
#SBATCH --array=0-1
#SBATCH --output=/results/vep_clinical/logs/build-%A_%a.log
set -euo pipefail

ASMS=(GRCh37 GRCh38)
ASM=${ASMS[$SLURM_ARRAY_TASK_ID]}
VC=/mnt/data/vcfcache/.venv/bin/vcfcache
SOURCE=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_${ASM}_joint_af001
DEST=/results/vep_clinical/databases/$ASM
NAME=vep_clinical_${ASM,,}_af001
RECIPE=/results/reproduction/final_clinical_campaign/vep_clinical_${ASM}.yaml
PARAM_SOURCE=$SOURCE/cache/vep115.2_everything/params.snapshot.yaml
WORK=/mnt/data/tmp/vep-clinical-build-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}
mkdir -p "$WORK" "$(dirname "$DEST")" /results/vep_clinical/logs
trap 'rm -rf "$WORK"' EXIT
sed -e "s#^temp_dir:.*#temp_dir: \"$WORK\"#" "$PARAM_SOURCE" > "$WORK/params.yaml"

for path in "$SOURCE/blueprint/vcfcache.bcf" "$SOURCE/blueprint/vcfcache.bcf.csi" \
            "$RECIPE" "$PARAM_SOURCE"; do
  [[ -s "$path" ]] || { echo "missing $path" >&2; exit 1; }
done

if [[ ! -s "$DEST/blueprint/vcfcache.bcf" ]]; then
  "$VC" blueprint-init --vcf "$SOURCE/blueprint/vcfcache.bcf" --output "$DEST"
fi

CACHE=$DEST/cache/$NAME
REUSE=false
if [[ -f "$CACHE/.vcfcache_complete" ]] && \
   cmp -s "$RECIPE" "$CACHE/annotation.yaml"; then
  REUSE=true
  echo "reusing completed cache with identical frozen recipe: $CACHE"
elif [[ -f "$CACHE/.vcfcache_complete" ]]; then
  echo "completed cache recipe differs; rebuilding $CACHE"
fi
if [[ "$REUSE" == false ]]; then
  /usr/bin/time -v -o "/results/vep_clinical/logs/build-${SLURM_ARRAY_JOB_ID}-${ASM}.time" \
    "$VC" cache-build --db "$DEST" --name "$NAME" \
      --anno-config "$RECIPE" --params "$WORK/params.yaml" --force
fi

BCF=$CACHE/vcfcache_annotated.bcf
[[ -s "$BCF" && -s "$BCF.csi" ]] || { echo "cache output missing" >&2; exit 1; }
HEADER=$(bcftools view -h "$BCF")
for field in CSQ CLIN_SIG gnomADe_AF; do
  grep -q "$field" <<< "$HEADER" || { echo "cache has no $field" >&2; exit 1; }
done
sha256sum "$SOURCE/blueprint/vcfcache.bcf" "$RECIPE" "$CACHE/annotation.yaml" \
  "$CACHE/params.snapshot.yaml" "$BCF" "$BCF.csi" \
  > "/results/vep_clinical/logs/build-${SLURM_ARRAY_JOB_ID}-${ASM}.sha256"
bcftools index -n "$BCF" \
  > "/results/vep_clinical/logs/build-${SLURM_ARRAY_JOB_ID}-${ASM}.records"
echo "completed $ASM $CACHE"
