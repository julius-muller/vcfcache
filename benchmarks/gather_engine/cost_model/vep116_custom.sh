#!/bin/bash
# The vep-rs mission tables (veprs_mission_stack.sh / _join.sh) as VCFs for Ensembl VEP's --custom.
# Ensembl VEP's dbNSFP plugin annotates only SNVs within transcripts, so the generic tables that
# vep-rs reads through its dbNSFP plugin (ClinVar, gnomAD exomes, dbSNP, PhyloP, mock tracks) go to
# Ensembl VEP as --custom VCFs with exact allele matching; same rows, same values.
# Writes /results/vep116/mission/custom/<name>.vcf.gz (+ .tbi) and touches /results/vep116/mission/CUSTOM.
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
ls /results /results/veprs05 /results/veprs05/mission/tables > /dev/null 2>&1
T=/results/veprs05/mission/tables; O=/results/vep116/mission/custom; mkdir -p $O
W=$(mktemp -d /mnt/data/tmp/vep116-custom-XXXX); trap 'rm -rf $W' EXIT; cd $W
conv() {  # table name
  local t=$1 n=$2 s=$(date +%s)
  local cols=$(zcat $T/$t.tsv.gz | head -1 | cut -f5-)
  { echo "##fileformat=VCFv4.2"
    for c in $cols; do echo "##INFO=<ID=$c,Number=.,Type=String,Description=\"$c (from $t.tsv.gz)\">"; done
    for c in $(seq 1 22) X Y M; do echo "##contig=<ID=chr$c>"; done
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    zcat $T/$t.tsv.gz | awk -F'\t' -v cols="$cols" 'BEGIN{OFS="\t"; n=split(cols,h," ")} NR>1 {
      s=""; for(i=1;i<=n;i++){ v=$(i+4); if(v=="."||v=="") continue; gsub(/[;= ]/,"_",v); s=s (s==""?"":";") h[i] "=" v }
      if(s=="") s="."; print $1,$2,".",$3,$4,".",".",s }'
  } | bgzip -@3 -l1 > $W/$n.vcf.gz && tabix -p vcf $W/$n.vcf.gz && mv $W/$n.vcf.gz $W/$n.vcf.gz.tbi $O/ \
    && echo "[$(date -Is)] $n: $(( $(date +%s) - s ))s, $(du -h $O/$n.vcf.gz | cut -f1)" || { echo "FAILED $n"; return 1; }
}
fail=0
for p in "clinvar clinvar" "gnomad gnomad" "dbsnp dbsnp" "phylop phylop" "trackA trackA" "trackB trackB"; do
  ( conv $p ) & done
for j in $(jobs -p); do wait $j || fail=1; done
[ $fail = 0 ] || { echo "FAILED"; exit 1; }
ls -la $O; touch $O/../CUSTOM; echo DONE
