#!/bin/bash
# Ensembl VEP 116.3 with the eleven sources of the vep-rs mission stack (vep116_stack.sh) on chr22 of
# one KPGP genome: does every source fill its fields, how does the fill compare with vep-rs 0.5.0
# (wgs_mission_veprs.sh) on the same input, and how fast is each (8 forks)?
# Needs vep116_custom.sh. Usage: probe_vep116.sh [sample]
set -uo pipefail
SAMPLE=${1:-KPGP-00128}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
ls /results /results/veprs05 /results/veprs05/mission /results/vep116/mission > /dev/null 2>&1
test -e /results/vep116/mission/CUSTOM || { echo "custom VCFs missing"; exit 1; }
VR=/results/veprs05/dist/vep-0.5.0-x86_64-unknown-linux-gnu/vep-0.5.0-x86_64-unknown-linux-gnu/vep
JC=/results/veprs/json_cache/GRCh38_native_20261023
T=/results/veprs05/mission/tables; CU=/results/vep116/mission/custom
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
OUT=/results/fastvep05/cost_model/probe_vep116; mkdir -p $OUT
W=$(mktemp -d /mnt/data/tmp/pv116-XXXX); trap 'rm -rf $W' EXIT; cd $W
VDIR=$W/vepdir; mkdir -p $VDIR; FORKS=8
. /results/fastvep05/cost_model/vep116_stack.sh
# vep-rs: its own mission stack (generic tables through its dbNSFP plugin)
D=$VD; Gr=""; for f in AF AC AN nhomalt AF_afr AF_amr AF_asj AF_eas AF_fin AF_mid AF_nfe AF_remaining AF_sas; do Gr="$Gr,gnomAD_$f"; done
PLR="--plugin dbNSFP,$T/clinvar.tsv.gz,CV_CLNSIG,CV_REVSTAT,CV_CLNDN --plugin dbNSFP,$T/gnomad.tsv.gz$Gr --plugin dbNSFP,$T/dbsnp.tsv.gz,dbSNP_RS"
PLR="$PLR --plugin REVEL,file=$D/new_tabbed_revel_grch38.tsv.gz --plugin AlphaMissense,$D/AlphaMissense_hg38.tsv.gz"
PLR="$PLR --plugin SpliceAI,snv=$T/spliceai_snv.vcf.gz,indel=$T/spliceai_indel.vcf.gz --plugin dbNSFP,$T/phylop.tsv.gz,PhyloP"
PLR="$PLR --plugin dbNSFP,$T/trackA.tsv.gz,TrackA_P,TrackA_Q --plugin dbNSFP,$T/trackB.tsv.gz,TrackB_P,TrackB_Q"
PLR="$PLR --plugin pLI,$D/pLI_values.txt --plugin LoFtool,$D/LoFtool_scores.txt --plugin GWAS,file=$D/gwas-catalog-download-associations-alt-full.tsv"
bcftools view -r chr22 $KP/$SAMPLE.GRCh38.small_variants.vcf.gz -Ov -o in.vcf
echo "node $(hostname); input $(grep -vc '^#' in.vcf) variants; $($VEP --help 2>&1 | grep -m1 ensembl-vep)" | tee $OUT/info.txt
# first run creates the GWAS plugin's processed file in --dir; the second is timed
for run in prime timed; do
  s=$(date +%s.%N); $VEP $VFLAGS -i in.vcf -o vep116.vcf --force_overwrite > vep116.log 2>&1; r=$?
  echo "vep116 $run $(echo "$(date +%s.%N) - $s" | bc) s, exit $r" | tee -a $OUT/info.txt
done
grep -iE "error|warn|fail" vep116.log | sort | uniq -c | head -20 | tee -a $OUT/info.txt
s=$(date +%s.%N); $VR --assembly GRCh38 --json_cache $JC --fork 8 --vcf --no_stats --format vcf -i in.vcf -o veprs.vcf \
  --force_overwrite --hgvs --fasta $FA $PLR > veprs.log 2>&1; r=$?
echo "veprs $(echo "$(date +%s.%N) - $s" | bc) s, exit $r" | tee -a $OUT/info.txt
python3 - vep116.vcf veprs.vcf <<'P' | tee $OUT/fill.tsv
import sys
def fill(path):
    fmt=None; filled={}; n=0
    for l in open(path):
        if l.startswith("##INFO=<ID=CSQ"): fmt=l.split("Format: ")[1].rstrip('">\n').split("|")
        if l.startswith("#"): continue
        n+=1; info=l.split("\t")[7]
        csq=[x for x in info.split(";") if x.startswith("CSQ=")]
        if not csq: continue
        got=set()
        for e in csq[0][4:].split(","):
            for k,v in zip(fmt,e.split("|")):
                if v: got.add(k)
        for k in got: filled[k]=filled.get(k,0)+1
    return n, fmt or [], filled
na, fa, a = fill(sys.argv[1]); nb, fb, b = fill(sys.argv[2])
print(f"field\tvep116 ({na})\tveprs ({nb})")
for k in fa + [x for x in fb if x not in fa]:
    print(f"{k}\t{a.get(k, '-' if k not in fa else 0)}\t{b.get(k, '-' if k not in fb else 0)}")
P
echo DONE
