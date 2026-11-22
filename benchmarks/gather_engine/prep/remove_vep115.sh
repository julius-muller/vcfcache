#!/bin/bash
# Remove the VEP 115 installation once 116.3 is verified on this node.
set -euo pipefail
A=/mnt/data/apps/ensembl-vep
test -e $A/116/READY.GRCh38 && test -e $A/116/READY.GRCh37
for asm in GRCh38 GRCh37; do
  for f in $A/115/cachedir/homo_sapiens/115_$asm/*.fa.gz*; do
    cmp "$f" "$A/116/cachedir/homo_sapiens/116_$asm/$(basename "$f")"
  done
done
test -d $A/115 && rm -rf $A/115
df -h /mnt/data | tail -1; ls $A
echo "$(hostname) removed VEP 115 at $(date -Is)" >> /results/apps-staging/logs/vep115-removed.txt
