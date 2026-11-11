#!/bin/bash
# Annotate one part of the blueprint with VEP 116.3, same flags as the cache recipe.
set -euo pipefail
ASM=$1; K=$2; N=$3
D=/mnt/data/apps/ensembl-vep/116
test -e $D/READY.$ASM
if [ "$ASM" = GRCh38 ]; then BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
else BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh37_joint_af001/blueprint/vcfcache.bcf; fi
test -s $BP
OUT=/results/vep116/build-$ASM; mkdir -p $OUT/aux-$K
# The container sees /mnt/data only, so VEP writes its side files node-locally.
AUX=$(mktemp -d /mnt/data/tmp/vep116-aux-XXXX); trap 'cp $AUX/* $OUT/aux-$K/ 2>/dev/null || true' EXIT
python3 /results/apps-staging/prep/parts.py $BP $N $K > $OUT/part-$K.bed
in=$(bcftools view -H -R $OUT/part-$K.bed --regions-overlap 0 $BP | wc -l)
echo "$(hostname) part $K/$N records $in contigs $(wc -l < $OUT/part-$K.bed)"
start=$(date +%s)
bcftools view -R $OUT/part-$K.bed --regions-overlap 0 $BP | \
  apptainer exec -B /mnt/data:/mnt/data $D/vep.sif vep --offline --buffer_size 100000 --fork 10 --cache \
    --dir_cache $D/cachedir --format vcf --vcf --force_overwrite -i stdin -o stdout -a $ASM \
    --hgvsg --everything --stats_file $AUX/vep_stats.html --warning_file $AUX/vep_warnings.txt \
    2> $AUX/vep_stderr.txt | \
  bcftools view -o $OUT/part-$K.tmp.bcf -Ob -W --threads 4
cp $AUX/* $OUT/aux-$K/ || true; rm -rf $AUX
out=$(bcftools index -n $OUT/part-$K.tmp.bcf)
echo "in=$in out=$out seconds=$(( $(date +%s) - start ))"
test "$out" -gt 0
mv $OUT/part-$K.tmp.bcf $OUT/part-$K.bcf; mv $OUT/part-$K.tmp.bcf.csi $OUT/part-$K.bcf.csi
echo "{\"part\": $K, \"node\": \"$(hostname)\", \"records\": $out, \"seconds\": $(( $(date +%s) - start ))}" > $OUT/part-$K.json
