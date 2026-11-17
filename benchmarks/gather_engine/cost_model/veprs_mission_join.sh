#!/bin/bash
# Join the per-chromosome parts of veprs_mission_stack.sh into genome-wide,
# tabix-indexed plugin tables, and build the ClinVar table from the same ClinVar
# release the fastVEP stack used. Usage: veprs_mission_join.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
ls /results /results/veprs05 > /dev/null 2>&1
OUT=/results/veprs05/mission; P=$OUT/parts; T=$OUT/tables; mkdir -p $T
CH="$(seq 1 22) X Y"
for c in $CH; do test -e $OUT/done/chr$c || { echo "chr$c missing"; exit 1; }; done
log() { echo "[$(date -Is)] $*"; }
table() {  # name parts header-columns...
  local n=$1 parts=$2; shift 2
  { (IFS=$'\t'; echo "#chr	pos(1-based)	ref	alt	$*"); for c in $CH; do [ -e $P/${parts}_chr$c.gz ] && zcat $P/${parts}_chr$c.gz; done; } | bgzip -@8 > $T/$n.tsv.gz && tabix -f -s1 -b2 -e2 -S1 $T/$n.tsv.gz
  log "$n: $(du -h $T/$n.tsv.gz | cut -f1)"
}
G="gnomAD_AF gnomAD_AC gnomAD_AN gnomAD_nhomalt"; for p in afr amr asj eas fin mid nfe remaining sas; do G="$G gnomAD_AF_$p"; done
table gnomad gnomad $G & table phylop phylop PhyloP & table dbsnp dbsnp dbSNP_RS & wait
table trackA track TrackA_P TrackA_Q & table trackB track TrackB_P TrackB_Q & wait
{ echo "##fileformat=VCFv4.2"; echo '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI max delta score from gnomAD v4.1: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">'
  for c in $CH; do echo "##contig=<ID=chr$c>"; done; echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  for c in $CH; do [ -e $P/spliceai_chr$c.gz ] && zcat $P/spliceai_chr$c.gz; done; } | bgzip -@8 > $T/spliceai.vcf.gz
bcftools view -v snps $T/spliceai.vcf.gz -Oz -W=tbi -o $T/spliceai_snv.vcf.gz && bcftools view -v indels $T/spliceai.vcf.gz -Oz -W=tbi -o $T/spliceai_indel.vcf.gz && rm $T/spliceai.vcf.gz
log "spliceai: $(du -ch $T/spliceai_*.vcf.gz | tail -1 | cut -f1)"
{ echo -e "#chr\tpos(1-based)\tref\talt\tCV_CLNSIG\tCV_REVSTAT\tCV_CLNDN"
  bcftools query -f 'chr%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\n' /results/fastvep05/acmg/sources/clinvar.vcf.gz | awk -F'\t' '$4!="."' ; } | bgzip -@8 > $T/clinvar.tsv.gz && tabix -f -s1 -b2 -e2 -S1 $T/clinvar.tsv.gz
log "clinvar: $(du -h $T/clinvar.tsv.gz | cut -f1)"
touch $OUT/TABLES; ls -la $T
