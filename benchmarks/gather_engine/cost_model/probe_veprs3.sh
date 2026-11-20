#!/bin/bash
# do the 1000G table and the real dbNSFP file fill values for chr-named input in vep-rs?
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
ls /results /results/sources /results/veprs05 > /dev/null 2>&1
V=/results/veprs05/dist/vep-0.5.0-x86_64-unknown-linux-gnu/vep-0.5.0-x86_64-unknown-linux-gnu/vep
T=/results/veprs/json_cache/GRCh38_native_20261023
W=$(mktemp -d /mnt/data/tmp/vp3-XXXX); trap 'rm -rf $W' EXIT; cd $W
DB=/results/sources/dbnsfp/dbNSFP5.4a_grch38.gz
bcftools view -H -r chr22:20000000-21000000 /results/sources/onekg/chr22.sites.vcf.gz | awk 'NR%20==0' | head -1500 | cut -f1-5 > b
tabix $DB 22:20000000-21000000 | awk 'NR%200==0' | head -500 | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,".",$3,$4}' >> b
sort -k2,2n b | awk '!s[$2" "$4" "$5]++' | awk '$5 ~ /^[ACGT]+$/' > b2
{ echo "##fileformat=VCFv4.2"; echo "##contig=<ID=chr22>"; echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; awk 'BEGIN{OFS="\t"}{print $1,$2,".",$4,$5,".",".","."}' b2; } > v.vcf
echo "input $(grep -vc '^#' v.vcf)"
{ echo -e "#chr\tpos(1-based)\tref\talt\tKG_AF\tKG_AF_EAS\tKG_AF_AMR\tKG_AF_EUR\tKG_AF_AFR\tKG_AF_SAS"
  bcftools norm -m- -r chr22:20000000-21000000 /results/sources/onekg/chr22.sites.vcf.gz -Ou 2>/dev/null | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AF_EAS\t%AF_AMR\t%AF_EUR\t%AF_AFR\t%AF_SAS\n'; } | bgzip > kg.tsv.gz; tabix -s1 -b2 -e2 -S1 kg.tsv.gz
$V --assembly GRCh38 --json_cache $T --fork 4 --vcf --no_stats --format vcf -i v.vcf -o o.vcf --force_overwrite \
  --plugin dbNSFP,kg.tsv.gz,KG_AF,KG_AF_EAS,KG_AF_AMR,KG_AF_EUR,KG_AF_AFR,KG_AF_SAS \
  --plugin dbNSFP,$DB,SIFT4G_score,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred \
  --plugin dbNSFP,$DB,GERP_92_mammals,phyloP100way_vertebrate > log 2>&1; echo "exit $?"; grep -iE "error|fail" log | head
python3 - o.vcf <<'P'
import sys
fmt=None; filled={}; n=0
for l in open(sys.argv[1]):
    if l.startswith("##INFO=<ID=CSQ"): fmt=l.split("Format: ")[1].rstrip('">\n').split("|")
    if l.startswith("#"): continue
    n+=1; info=l.split("\t")[7]
    csq=[x for x in info.split(";") if x.startswith("CSQ=")][0][4:]
    got=set()
    for e in csq.split(","):
        for k,v in zip(fmt,e.split("|")):
            if k.startswith("dbNSFP") and v: got.add(k)
    for k in got: filled[k]=filled.get(k,0)+1
print("records", n); [print(k, v) for k,v in sorted(filled.items())]
P
echo DONE
