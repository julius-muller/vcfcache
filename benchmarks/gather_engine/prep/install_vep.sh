#!/bin/bash
# Install VEP 116.3 (image, offline cache, FASTA) node-locally, beside the 115 installation.
set -euo pipefail
ASM=$1
ST=/results/apps-staging/ensembl-vep-116
D=/mnt/data/apps/ensembl-vep/116
OLD=/mnt/data/apps/ensembl-vep/115/cachedir/homo_sapiens/115_$ASM
NEW=$D/cachedir/homo_sapiens/116_$ASM
test -e $ST/DOWNLOADED.$ASM
# Scratch of this session's own cancelled benchmark jobs.
find /mnt/data/tmp -maxdepth 1 -user "$USER" \( -name 'cohort-*' -o -name 'gather-*' \) -exec rm -rf {} + || true
df -h /mnt/data | tail -1
mkdir -p $D/cachedir
if [ ! -e $D/READY.$ASM ]; then
  # A completed extraction is recognised by the FASTA index copied after it.
  if ! ls $NEW/*.fa.gz.gzi >/dev/null 2>&1; then
    rm -rf $NEW
    tar xzf $ST/homo_sapiens_vep_116_$ASM.tar.gz -C $D/cachedir
    test -s $NEW/info.txt
    for f in $OLD/*.fa.gz $OLD/*.fa.gz.fai $OLD/*.fa.gz.gzi; do cp --reflink=auto $f $NEW/; done
  fi
  # No -p: the shared mount carries ACLs the local disk cannot take.
  cp $ST/vep.sif $D/vep.sif.tmp && chmod 755 $D/vep.sif.tmp && mv $D/vep.sif.tmp $D/vep.sif
fi
sha256sum -c <(sed "s#  vep.sif#  $D/vep.sif#" $ST/vep.sif.sha256)
V="apptainer exec -B /mnt/data:/mnt/data $D/vep.sif vep"
$V --help | grep -i "ensembl-vep"
# Functional check: offline --everything on 300 blueprint records.
if [ "$ASM" = GRCh38 ]; then BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
else BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh37_joint_af001/blueprint/vcfcache.bcf; fi
W=$(mktemp -d /mnt/data/tmp/vep116-check-XXXX)
set +o pipefail; { bcftools view -h $BP; bcftools view -H $BP | head -n 300; } > $W/in.vcf; set -o pipefail
$V --offline --fork 4 --cache --dir_cache $D/cachedir --format vcf --vcf --force_overwrite -i $W/in.vcf -o $W/out.vcf -a $ASM --hgvsg --everything --no_stats 2> $W/err.txt || { tail -20 $W/err.txt; exit 1; }
n=$(grep -v "^#" $W/out.vcf | grep -c "CSQ="); echo "annotated records: $n"; test "$n" -gt 0
grep "^##VEP=" $W/out.vcf | cut -c1-300
rm -rf $W
cat > $D/node_ready.$ASM.json <<J
{"hostname": "$(hostname)", "assembly": "$ASM", "vep": "$($V --help | awk '/ensembl-vep/ {print $NF}')", "sif_sha256": "$(sha256sum $D/vep.sif | cut -d' ' -f1)", "cache": "$NEW", "fasta_from": "$OLD", "date": "$(date -Is)"}
J
touch $D/READY.$ASM; cat $D/node_ready.$ASM.json; du -sh $NEW; df -h /mnt/data | tail -1
