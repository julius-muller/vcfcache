#!/bin/bash
# One shared set of SnpSift sources for wgs_mission_snpeff.sh, so that no node builds its own index.
# Run on a node whose /mnt/data/tmp/snpeff-mission/src holds the decompressed sources with their
# .sidx (snpsift_memory_probe.sh built all but dbSNP's there). dbSNP's .sidx did not fit a 26 GB heap
# with SnpSift's default block size (16384); this tries a block size of 1 MB. If that builds and
# annotation through it stays within 6 GB, dbSNP is read as plain VCF like the others; otherwise
# through its tabix index (same output, 15-20 times slower, snpsift_index_probe.sh).
# Writes /results/snpeff/mission/src and touches READY there.
# Usage: snpeff_sources.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
M=/results/snpeff/mission; S=$M/src; mkdir -p $S
R=/mnt/data/tmp/snpeff-mission; SRC=$R/src
W=$(mktemp -d /mnt/data/tmp/ses-XXXX); trap 'rm -rf $W' EXIT; cd $W
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname)"
ERR=$W/err; CU=$SRC; T=$SRC; . /results/fastvep05/cost_model/snpeff_stack.sh
for s in clinvar gnomad phylop spliceai_snv spliceai_indel trackA trackB; do
  test -s $SRC/$s.vcf.sidx || { echo "missing $SRC/$s.vcf.sidx"; exit 1; }; done
SIFT="$JAVA -Xmx26g -jar $SE/SnpSift.jar annotate -noLog -c $SE/snpEff/snpEff.config -noId"
bcftools view -h $SRC/clinvar.vcf > one.vcf; bcftools view -H $SRC/clinvar.vcf | head -1 >> one.vcf
bcftools view -r chr22 $KP/KPGP-00128.GRCh38.small_variants.vcf.gz -Ov | eval "$SNPEFF" > se.vcf
m() {  # step command: prints wall and peak RSS (GB); fails if the command fails
  /usr/bin/time -v -o tv bash -c "$2" 2>> $ERR; local rc=$?
  printf "%s\t%s\t%s\t%s\n" "$1" "$(awk -F': ' '/Elapsed/ {print $2}' tv)" \
    "$(awk -F': ' '/Maximum resident/ {printf "%.2f", $2/1048576}' tv)" "$rc" | tee -a $M/memory_probe/probe.tsv
  return $rc
}
rm -f $SRC/dbsnp.vcf.sidx
MODE=tabix
if m build_dbsnp_26g_block1m "$SIFT -maxBlockSize 1048576 $SRC/dbsnp.vcf one.vcf > /dev/null" \
   && m annotate_dbsnp_26g_block1m "$SIFT $SRC/dbsnp.vcf se.vcf > plain.vcf"; then
  rss=$(tail -n 1 $M/memory_probe/probe.tsv | cut -f3)
  $SIFT /results/vep116/mission/custom/dbsnp.vcf.gz se.vcf > tabix.vcf 2>> $ERR
  same=$(cmp -s <(grep -v '^##SnpSift' plain.vcf) <(grep -v '^##SnpSift' tabix.vcf) && echo yes || echo no)
  log "dbSNP plain: annotation peak ${rss} GB, output same as tabix: $same"
  [ "$same" = yes ] && awk -v r=$rss 'BEGIN{exit !(r <= 6)}' && MODE=plain
else tail -n 3 $ERR; fi
log "dbSNP read as $MODE"
rm -f $S/READY
for s in clinvar gnomad phylop spliceai_snv spliceai_indel trackA trackB; do cp -p $SRC/$s.vcf $SRC/$s.vcf.sidx $S/ || exit 1; done
if [ $MODE = plain ]; then cp -p $SRC/dbsnp.vcf $SRC/dbsnp.vcf.sidx $S/ || exit 1
else cp -p /results/vep116/mission/custom/dbsnp.vcf.gz /results/vep116/mission/custom/dbsnp.vcf.gz.tbi $S/ || exit 1; fi
echo $MODE > $S/DBSNP_MODE; touch $S/READY
log "shared sources: $(du -sh $S | cut -f1)"; echo DONE
