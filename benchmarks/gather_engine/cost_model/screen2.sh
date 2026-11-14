#!/bin/bash
# Stage 1 screen, part 2.
# A: does the fixed cost of the ACMG and dense stacks grow with the span of
#    genome the input covers (source blocks) or is it start-up?
# B: real SpliceAI2 precomputed scores (84-column table) as a fastVEP custom
#    source on a few megabases of chr21: all score columns, and five of them.
# Usage: screen2.sh [START END]   (chr21 region for part B)
set -uo pipefail
LO=${1:-32000000}; HI=${2:-37000000}
OUT=/results/fastvep05/cost_model/screen2; mkdir -p $OUT
W=/mnt/data/tmp/screen2; rm -rf $W; mkdir -p $W; cd $W
SRC=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp/KPGP-00128.GRCh38.small_variants.vcf.gz
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
FV=/mnt/data/apps/fastvep/0.5.0
SA=/results/fastvep05/acmg/sa
S2=/mnt/data/tmp/spliceai2/GRCh38
BASE="--transcript-cache $FV/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache --fasta $FV/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa --no-progress --hgvs"
: > $OUT/runs.tsv
run() {  # part stack input flags...
  local part=$1 s=$2 i=$3; shift 3
  local s0=$(date +%s.%N)
  /usr/bin/time -f "%U\t%S\t%M" -o t.txt bash -c "bcftools view $i -Ov | env RAYON_NUM_THREADS=8 $FV/bin/fastvep annotate --input - --output - --output-format vcf $BASE $* 2> err.txt | bcftools view -Ob -o d.bcf -W --threads 8" || tail -2 err.txt
  local e=$(echo "$(date +%s.%N) - $s0" | bc)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $part $s $i "$(bcftools index -n $i)" "$e" "$(tail -1 t.txt)" "$(bcftools index -n d.bcf)" "$(stat -c %s d.bcf)" | tee -a $OUT/runs.tsv
}
thin() {  # in out k
  { bcftools view -h $1; bcftools view -H $1 | awk -v k=$3 'NR % k == 0'; } | bcftools view -Ob1 -W -o $2
}
# ---- A
bcftools view $SRC -r chr1 -Ob1 -W --threads 8 -o c1.bcf
thin c1.bcf wide4.bcf 4
bcftools view c1.bcf -r chr1:1-62000000 -Ob1 -W -o narrow.bcf
for rep in 1 2; do
  for s in acmg:"--acmg --sa-dir $SA/stack" dense:"--sa-dir $SA/dense" plain:""; do
    run A ${s%%:*} wide4.bcf ${s#*:}; run A ${s%%:*} narrow.bcf ${s#*:}
  done
done
# ---- B
names=$(zcat $S2/snv/chr21.tsv.gz | head -1 | cut -f6-)
{
  echo "##fileformat=VCFv4.2"; echo "##contig=<ID=chr21,length=46709983>"
  echo '##INFO=<ID=strand,Number=1,Type=String,Description="strand">'
  for n in $names; do echo "##INFO=<ID=$n,Number=1,Type=Float,Description=\"SpliceAI2 $n\">"; done
  printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
} > hdr.txt
conv() {  # kind
  zcat $S2/$1/chr21.tsv.gz | awk -F'\t' -v lo=$LO -v hi=$HI 'NR == 1 { for (i = 5; i <= NF; i++) name[i] = $i; next }
    $2 > hi { exit } $2 >= lo { s = "strand=" $5; for (i = 6; i <= NF; i++) if ($i != "" && $i != ".") s = s ";" name[i] "=" $i
      print $1 "\t" $2 "\t.\t" $3 "\t" $4 "\t.\t.\t" s }'
}
s0=$(date +%s)
{ cat hdr.txt; { conv snv; conv indel; } | sort -t$'\t' -k2,2n -S 4G --stable; } | bgzip -@8 > sai2.vcf.gz
echo "converted $(( $(date +%s) - s0 ))s size $(stat -c %s sai2.vcf.gz) records $(zcat sai2.vcf.gz | grep -vc '^#')" | tee $OUT/build.txt
few="spliceai2_summary_score,donor_gain_delta_score_0,donor_loss_delta_score_0,acceptor_gain_delta_score_0,acceptor_loss_delta_score_0"
mkdir -p sa_all sa_few
/usr/bin/time -f "build all: %e s, %M KB" $FV/bin/fastvep sa-build --source custom_vcf --name SpliceAI2 -i sai2.vcf.gz -o sa_all/spliceai2 --chunk-bits 16 2>&1 | tail -3 | tee -a $OUT/build.txt
/usr/bin/time -f "build few: %e s, %M KB" $FV/bin/fastvep sa-build --source custom_vcf --name SpliceAI2 --info-fields $few -i sai2.vcf.gz -o sa_few/spliceai2 --chunk-bits 16 2>&1 | tail -3 | tee -a $OUT/build.txt
ls -la sa_all sa_few | tee -a $OUT/build.txt
R=chr21:$LO-$HI
bcftools view $SRC -r $R -Ob1 -W -o smp.bcf
bcftools view $BP -r $R -Ob1 -W -o bp1.bcf
for k in 2 4 8; do thin bp1.bcf bp$k.bcf $k; done
for rep in 1 2; do
  for i in bp1 bp2 bp4 bp8 smp; do
    run B plain $i.bcf ""; run B sai2_few $i.bcf --sa-dir sa_few; run B sai2_all $i.bcf --sa-dir sa_all
  done
done
bcftools view -H d.bcf | awk -F'\t' '{print length($8)}' | sort -n | awk '{a[NR]=$1} END{print "INFO bytes per record with all columns: median", a[int(NR/2)], "max", a[NR]}' | tee -a $OUT/build.txt
bcftools view -H d.bcf | grep -c "SpliceAI2" | sed 's/^/records carrying SpliceAI2: /' | tee -a $OUT/build.txt
cd /; rm -rf $W; echo DONE
