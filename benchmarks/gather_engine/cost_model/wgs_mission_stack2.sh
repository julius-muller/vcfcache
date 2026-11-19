#!/bin/bash
# Five further well-known sources for the mission stack, per chromosome, at
# fastVEP's advised width (65 kb, sa16) and narrow (sa12, misses only):
#   gnomAD v4.1 genomes (custom source gnomAD_genomes: fastVEP keeps one answer
#   per allele under the built-in `gnomad` key, which the exomes already hold),
#   1000 Genomes 30x panel (onekg), GERP 91 mammals (gerp), dbNSFP 5.4a (dbnsfp)
#   and, once, ENCODE cCREs (custom BED).
# Usage: wgs_mission_stack2.sh <chromosome without chr>
set -uo pipefail
C=$1
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
FV=/mnt/data/apps/fastvep/0.5.0/bin/fastvep
SRC=/results/sources; TOOLS=/results/fastvep05/tools
ls /results /results/fastvep05 $SRC > /dev/null 2>&1   # clears a stale "no such directory" on some nodes
# NB=<bits>: narrow copies at another width (default 12) into /results/fastvep05/mission2<NB>/sa12.
# gnomAD genomes only with GG=1: as 20 Mb custom pieces it is not a valid fastVEP source (see AGENTS.md).
NB=${NB:-12}; GG=${GG:-0}
OUT=/results/fastvep05/mission2${NB#12}; mkdir -p $OUT/sa16 $OUT/sa12 $OUT/done
[ -e $OUT/done/chr$C ] && exit 0
W=$(mktemp -d /mnt/data/tmp/mission2-$C-XXXX); trap 'rm -rf $W' EXIT; cd $W; mkdir 16 12
log() { echo "[$(date -Is)] chr$C $*"; }
b() { local s=$(date +%s); $FV sa-build --no-progress --assembly GRCh38 "$@" > $W/b.$BASHPID.log 2>&1 || { cat $W/b.$BASHPID.log; echo FAILED "$@"; return 1; }; log "sa-build $2 ${*: -3} $(( $(date +%s) - s ))s"; }
both() {  # source input name [extra...]: build at 16 and 12 bits, in parallel
  local src=$1 in=$2 name=$3; shift 3
  b --source $src -i $in -o $W/16/$name "$@" --chunk-bits 16 & local p=$!
  b --source $src -i $in -o $W/12/$name "$@" --chunk-bits $NB || { wait $p; return 1; }
  wait $p
}
fail=0
# gnomAD genomes: download, keep the fields fastVEP's gnomad source keeps, split into 20 Mb pieces
[ $GG != 1 ] || ( G=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr$C.vcf.bgz
  curl -fsSL --retry 10 --retry-delay 30 -o gg.vcf.bgz $G && curl -fsSL --retry 10 -o gg.vcf.bgz.tbi $G.tbi || { echo "FAILED gnomAD genomes download"; exit 1; }
  log "gnomAD genomes downloaded $(du -h gg.vcf.bgz | cut -f1)"
  F="AF,AC,AN,nhomalt"; for p in afr amr asj eas fin mid nfe remaining sas; do F="$F,AF_$p"; done
  K=$(echo $F | sed 's/[^,]*/INFO\/&/g')
  bcftools annotate -x "^$K" gg.vcf.bgz -Oz -o slim.vcf.gz -W=tbi --threads 4 && rm gg.vcf.bgz || exit 1
  len=$(bcftools view -h slim.vcf.gz | grep -o "ID=chr$C,length=[0-9]*" | cut -d= -f3)
  i=0; for ((s=1; s<=len; s+=20000000)); do
    i=$((i+1)); e=$((s+19999999))
    bcftools view -r chr$C:$s-$e slim.vcf.gz -Oz -o p.vcf.gz || exit 1
    [ $(bcftools view -H p.vcf.gz | head -1 | wc -l) = 0 ] && continue
    both custom_vcf p.vcf.gz gnomad_genomes_chr${C}_p$i --name gnomAD_genomes --info-fields $F || exit 1
  done ) & p1=$!
# 1000 Genomes 30x panel, sites only (no chrY panel)
( if [ $C != Y ]; then
    v=""; [ $C = X ] && v=".v2"
    K=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr$C.filtered.SNV_INDEL_SV_phased_panel$v.vcf.gz
    curl -fsSL --retry 10 --retry-delay 30 $K | bcftools view -G -Oz -o kg.vcf.gz || { echo "FAILED 1000 Genomes"; exit 1; }
    both onekg kg.vcf.gz onekg_chr$C
  fi ) & p2=$!
# GERP: bigWig -> bedGraph with chr names
( $TOOLS/bigWigToBedGraph -chrom=$C $SRC/gerp/gerp_91_mammals.GRCh38.bw gerp.raw.bg &&
  awk 'BEGIN{OFS="\t"} {print "chr"$1,$2,$3,$4}' gerp.raw.bg > gerp.bg && rm gerp.raw.bg &&
  both gerp gerp.bg gerp_chr$C ) & p3=$!
# dbNSFP 5.4a: this chromosome's rows with the header
( D=$SRC/dbnsfp/dbNSFP5.4a_grch38.gz
  { tabix -H $D; tabix $D $C; } | bgzip -@2 > dbn.tsv.gz && [ $(zcat dbn.tsv.gz | head -2 | wc -l) = 2 ] || { echo "no dbNSFP rows"; exit 0; }
  both dbnsfp dbn.tsv.gz dbnsfp_chr$C ) & p4=$!
if [ $C = 1 ]; then
  b --source custom_bed --name ENCODE_cCRE -i $SRC/encode/GRCh38-cCREs.bed -o $W/16/encode_ccre && cp $W/16/encode_ccre.osi $W/12/ || fail=1
fi
for p in $p1 $p2 $p3 $p4; do wait $p || fail=1; done
[ $GG = 1 ] && { ls $W/16/gnomad_genomes_chr${C}_p*.osa2 > /dev/null 2>&1 || { echo "missing gnomad_genomes"; fail=1; } }
for f in gerp $( [ $C != Y ] && echo onekg ); do
  [ -s $W/16/${f}_chr$C.osa2 ] && [ -s $W/12/${f}_chr$C.osa2 ] || { echo "missing $f"; fail=1; }
done
[ $(ls $W/16 | wc -l) = $(ls $W/12 | wc -l) ] || { echo "source sets differ"; fail=1; }
[ $fail = 0 ] || { echo "chr$C incomplete"; exit 1; }
for try in 1 2 3 4 5; do
  ls $(dirname $OUT) $OUT $OUT/sa16 $OUT/sa12 > /dev/null 2>&1; cp $W/16/* $OUT/sa16/ && cp $W/12/* $OUT/sa12/ && break; sleep 20
done
[ $(ls $OUT/sa16 | grep -c "_chr${C}[._]") -ge $(ls $W/16 | grep -c "_chr${C}[._]") ] || { echo "chr$C not saved"; exit 1; }
touch $OUT/done/chr$C
log "done: $(ls $W/16 | tr '\n' ' ')"
