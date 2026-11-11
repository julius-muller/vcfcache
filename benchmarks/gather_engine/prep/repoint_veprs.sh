#!/bin/bash
# Point the vep-rs clinical recipe at the FASTA under the 116 installation (same bytes as before).
set -euo pipefail
test "$(wc -l < /results/apps-staging/logs/vep115-removed.txt)" -ge 6
for f in /results/veprs05/caches/veprs05_clinical/cache/veprs05_clinical/annotation.yaml /results/veprs05/build_cache.sh; do
  cp -p $f $f.before-vep116
  sed -i 's#/mnt/data/apps/ensembl-vep/115/cachedir/homo_sapiens/115_GRCh38/#/mnt/data/apps/ensembl-vep/116/cachedir/homo_sapiens/116_GRCh38/#g' $f
  grep -c "ensembl-vep/116/cachedir/homo_sapiens/116_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz" $f
  ! grep -q "ensembl-vep/115" $f
done
