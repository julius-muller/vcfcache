#!/bin/bash
# Install fastVEP $V node-locally with rebuilt transcript caches.
set -euo pipefail
V=${1:?fastVEP version, e.g. 0.5.0}
ST=/results/apps-staging/fastvep
OLD=/mnt/data/apps/fastvep/0.3.0-publication/data
D=/mnt/data/apps/fastvep/$V
test -e $ST/DOWNLOADED.$V
mkdir -p $D/bin $D/data/GRCh38 $D/data/GRCh37
tar xzf $ST/fastVEP-v$V-x86_64-unknown-linux-gnu.tar.gz -C $D --strip-components=1
mv -f $D/fastvep $D/bin/fastvep
FV=$D/bin/fastvep; $FV --version
# Reference files are hard links to the 0.3.0 copies (same bytes, no extra space).
ln -f $OLD/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa $OLD/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.fai $D/data/GRCh38/
ln -f $OLD/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa $OLD/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai $OLD/GRCh37/Homo_sapiens.GRCh37.87.gff3 $D/data/GRCh37/
gzip -dc $ST/Homo_sapiens.GRCh38.116.gff3.gz > $D/data/GRCh38/Homo_sapiens.GRCh38.116.gff3
$FV cache --gff3 $D/data/GRCh38/Homo_sapiens.GRCh38.116.gff3 --fasta $D/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa -o $D/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache --no-progress
$FV cache --gff3 $D/data/GRCh37/Homo_sapiens.GRCh37.87.gff3 --fasta $D/data/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -o $D/data/GRCh37/Homo_sapiens.GRCh37.fastvep.cache --no-progress
W=$(mktemp -d /mnt/data/tmp/fv05-check-XXXX)
for ASM in GRCh38 GRCh37; do
  if [ $ASM = GRCh38 ]; then BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf; FA=$D/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
  else BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh37_joint_af001/blueprint/vcfcache.bcf; FA=$D/data/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa; fi
  set +o pipefail; { bcftools view -h $BP; bcftools view -H $BP | head -n 2000; } > $W/in.vcf; set -o pipefail
  $FV annotate --input $W/in.vcf --output $W/out.vcf --output-format vcf --transcript-cache $D/data/$ASM/Homo_sapiens.$ASM.fastvep.cache --fasta $FA --hgvs --no-progress --sa-dir /results/richness/sa/clinical_$ASM 2> $W/err.txt || { tail $W/err.txt; exit 1; }
  echo "$ASM annotated $(grep -v '^#' $W/out.vcf | grep -c 'CSQ=') with gnomAD $(grep -c 'FV_GNOMAD=' $W/out.vcf || true)"
  test "$(grep -v '^#' $W/out.vcf | grep -c 'CSQ=')" -gt 0
done
rm -rf $W
cat > $D/node_ready.json <<J
{"hostname": "$(hostname)", "fastvep_version": "$($FV --version)", "fastvep_sha256": "$(sha256sum $FV | cut -d' ' -f1)",
 "GRCh38": {"gff": "Homo_sapiens.GRCh38.116.gff3", "transcript_cache_sha256": "$(sha256sum $D/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache | cut -d' ' -f1)"},
 "GRCh37": {"gff": "Homo_sapiens.GRCh37.87.gff3", "transcript_cache_sha256": "$(sha256sum $D/data/GRCh37/Homo_sapiens.GRCh37.fastvep.cache | cut -d' ' -f1)"},
 "complete": true, "date": "$(date -Is)"}
J
cat $D/node_ready.json; du -sh $D
