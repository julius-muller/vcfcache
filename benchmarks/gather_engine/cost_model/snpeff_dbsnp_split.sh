#!/bin/bash
# dbSNP for SnpSift as five plain VCFs by chromosome group, each with its own .sidx. SnpSift could not
# build the .sidx of the whole file (1.1 billion records) in a 26 GB heap, with its default block size
# or with 1 MB blocks (snpeff_sources.sh), and through tabix it is 15-20 times slower. Each record of
# an input matches at most one group, so five SnpSift stages give the same output as one; the check
# below compares chr22 and 20 Mb of chr1, chr10 and chrX of KPGP-00128 against SnpSift through the tabix index of the whole file.
# Adds the groups to /results/snpeff/mission/src and sets DBSNP_MODE to split.
# Usage: snpeff_dbsnp_split.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
M=/results/snpeff/mission; S=$M/src; DB=/results/vep116/mission/custom/dbsnp.vcf.gz
test -e $S/READY || { echo "shared sources missing"; exit 1; }
# groups already built on this node (a check that failed) are kept
W=/mnt/data/tmp/snpeff-dbsnp; mkdir -p $W; cd $W
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname)"
ERR=$W/err; CU=$S; T=$S; . /results/fastvep05/cost_model/snpeff_stack.sh
$JAVA -version 2>&1 | head -1
SIFT="$JAVA -Xmx26g -jar $SE/SnpSift.jar annotate -noLog -c $SE/snpEff/snpEff.config -noId"
G1="chr1 chr2 chr3"; G2="chr4 chr5 chr6 chr7"; G3="chr8 chr9 chr10 chr11 chr12"; G4="chr13 chr14 chr15 chr16 chr17 chr18"
G5=$(tabix -l $DB | grep -v -x -F -f <(echo $G1 $G2 $G3 $G4 | tr ' ' '\n') | tr '\n' ' ')
[ -n "$G5" ] || { echo "FAILED: no contigs for group 5"; exit 1; }
bcftools view -h $DB > one.vcf; bcftools view -H $DB chr1 | head -1 >> one.vcf
for g in 1 2 3 4 5; do
  eval "regs=\$G$g"; s=$(date +%s)
  [ dbsnp_g$g.vcf.sidx -nt dbsnp_g$g.vcf ] && { log "group $g ($regs): kept"; continue; }
  bcftools view $DB -r $(echo $regs | tr ' ' ',') -Ov -o dbsnp_g$g.vcf || { echo "FAILED extract $g"; exit 1; }
  /usr/bin/time -v -o tv $SIFT dbsnp_g$g.vcf one.vcf > /dev/null 2>> $ERR || { echo "FAILED sidx $g"; tail -n 3 $ERR; exit 1; }
  log "group $g ($regs): $(du -h dbsnp_g$g.vcf | cut -f1), $(( $(date +%s) - s )) s, index peak $(awk -F': ' '/Maximum resident/ {printf "%.1f", $2/1048576}' tv) GB"
done
# same output as one dbSNP stage through tabix?
for c in chr22 chr1:1-20000000 chr10:1-20000000 chrX:1-20000000; do
  bcftools view -r $c $KP/KPGP-00128.GRCh38.small_variants.vcf.gz -Ov | eval "$SNPEFF" > se.vcf
  s=$(date +%s); $SIFT $DB se.vcf 2>> $ERR | grep -v '^#' > tabix.vcf; e1=$(( $(date +%s) - s ))
  s=$(date +%s)
  cat se.vcf | $SIFT dbsnp_g1.vcf - 2>> $ERR | $SIFT dbsnp_g2.vcf - 2>> $ERR | $SIFT dbsnp_g3.vcf - 2>> $ERR \
    | $SIFT dbsnp_g4.vcf - 2>> $ERR | $SIFT dbsnp_g5.vcf - 2>> $ERR | grep -v '^#' > split.vcf; e2=$(( $(date +%s) - s ))
  cmp -s tabix.vcf split.vcf || { echo "FAILED: $c differs"; diff tabix.vcf split.vcf | head -5; exit 1; }
  log "$c: $(wc -l < split.vcf) records identical; tabix ${e1} s, five plain groups ${e2} s"
done
rm -f $S/READY
for g in 1 2 3 4 5; do cp --preserve=timestamps dbsnp_g$g.vcf dbsnp_g$g.vcf.sidx $S/ || exit 1; done
# the copy can give a VCF and its index the same time stamp, and SnpSift rebuilds an index that is not newer
for g in 1 2 3 4 5; do touch $S/dbsnp_g$g.vcf.sidx; done
rm -f $S/dbsnp.vcf.gz $S/dbsnp.vcf.gz.tbi; echo split > $S/DBSNP_MODE; touch $S/READY
rm -rf $W; log "shared sources: $(du -sh $S | cut -f1)"; echo DONE
