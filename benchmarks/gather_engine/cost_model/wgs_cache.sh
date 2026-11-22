#!/bin/bash
# Annotate the whole gnomAD blueprint with fastVEP's ACMG stack to make the
# genome-wide vcfcache cache. Usage: wgs_cache.sh [blueprint tag, af001]
set -euo pipefail
TAG=${1:-af001}
D=/results/fastvep05/acmg; SA=$D/sa/wgs
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
for c in $(seq 1 22) X Y; do test -e $SA/chr$c.DONE; done
for f in clinvar.osa2 clinvar_protein.oga gnomad_genes.oga omim.oga repeatmasker.osi; do cp -n $D/sa/stack/$f $SA/; done
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_$TAG/blueprint/vcfcache.bcf
OLD=/results/fastvep05/caches/fastvep050_clinical_grch38_af001/cache/fastvep050_clinical_grch38_af001
NAME=fastvep050_acmg_grch38_$TAG; ROOT=/results/fastvep05/caches/$NAME; C=$ROOT/cache/$NAME
[ -e $ROOT.COMPLETE ] && exit 0
rm -rf $ROOT; mkdir -p $C $ROOT/blueprint; cp -r $(dirname $(dirname $OLD))/workflow $ROOT/
sed -e "s#--sa-dir /results/richness/sa/clinical_GRCh38#--sa-dir $SA --acmg#" $OLD/annotation.yaml > $C/annotation.yaml
grep -q -- "--sa-dir $SA --acmg" $C/annotation.yaml
cp $OLD/params.snapshot.yaml $C/
W=$(mktemp -d /mnt/data/tmp/wgs-cache-XXXX); trap 'rm -rf $W' EXIT
start=$(date +%s)
bcftools view $BP -Ov | env RAYON_NUM_THREADS=16 $FV annotate --input - --output - --output-format vcf \
  --transcript-cache $FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache --fasta $FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa \
  --hgvs --no-progress --sa-dir $SA --acmg 2> $W/err.txt | bcftools view -Ob -o $W/out.bcf -W --threads 8
in=$(bcftools index -n $BP); out=$(bcftools index -n $W/out.bcf); tail -5 $W/err.txt
test "$in" = "$out"
cp $W/out.bcf $C/vcfcache_annotated.bcf; cp $W/out.bcf.csi $C/vcfcache_annotated.bcf.csi
echo "{\"stack\": \"acmg\", \"sa_dir\": \"$SA\", \"records\": $out, \"seconds\": $(( $(date +%s) - start )), \"bytes\": $(stat -c %s $C/vcfcache_annotated.bcf)}" | tee $C/build_provenance.json
touch $ROOT.COMPLETE
