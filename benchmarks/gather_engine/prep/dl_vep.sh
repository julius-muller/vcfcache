#!/bin/bash
# Download the Ensembl VEP 116 offline cache for one assembly (and the 116.3 image) once, to shared storage.
set -euo pipefail
ASM=$1
ST=/results/apps-staging/ensembl-vep-116; mkdir -p $ST; cd $ST
U=https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache
T=homo_sapiens_vep_116_$ASM.tar.gz
curl -fsS -o CHECKSUMS.$ASM $U/CHECKSUMS
want=$(awk -v f=$T '$3==f {print $1" "$2}' CHECKSUMS.$ASM)
for try in 1 2 3 4 5 6; do
  curl -fS --retry 10 --retry-delay 30 -C - -o $T.part $U/$T && break || echo "retry $try"
done
have=$(sum $T.part | awk '{print $1" "$2}')
echo "checksum want=[$want] have=[$have]"
test "$want" = "$have"
mv $T.part $T
if [ "$ASM" = GRCh38 ] && [ ! -s vep.sif ]; then
  export APPTAINER_CACHEDIR=/mnt/data/tmp/apptainer-cache-$USER APPTAINER_TMPDIR=/mnt/data/tmp
  mkdir -p $APPTAINER_CACHEDIR
  apptainer pull --force vep.sif.part docker://ensemblorg/ensembl-vep:release_116.3
  mv vep.sif.part vep.sif
  apptainer exec vep.sif vep --help | grep -i "ensembl-vep" | tee vep.version
fi
sha256sum $T > $T.sha256
[ -s vep.sif ] && sha256sum vep.sif > vep.sif.sha256 || true
touch $ST/DOWNLOADED.$ASM
ls -la $ST
