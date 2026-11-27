#!/bin/bash
# Prep of the 40-genome campaign (cohort40.sh): the 20 KPGP and 20 SGDP evaluation genomes (GRCh38)
# of the external WGS manifest, a shared cohort cache from its 6 GRCh38 training genomes (3 KPGP,
# 3 SGDP), and the miss rates of every genome, including the 15 PGP genomes (GRCh37, miss rate only).
# Writes to /results/cohort40/prep:
#   samples.tsv              cohort, sample, role, superpopulation, assembly, path, records
#   sites/<sample>.bcf       biallelic sites of each GRCh38 genome (bcftools norm -m-, no INFO)
#   train.bcf                sites of the 6 training genomes
#   new46.bcf                sites of all 46 GRCh38 genomes beyond the AF >= 1% blueprint
#   extra.bcf                the part of new46 beyond the 23 KPGP genomes, whose sites the VEP 116.3
#                            and SnpEff missions annotated already (/results/vep116/mission/parts/new.bcf)
#   {bp,new,extra}-<k>.regions   parts by position: bp and new as in the VEP 116.3 mission
#   missrate.tsv             per genome: sites, beyond the blueprint, beyond blueprint and training set
# Usage: cohort40_prep.sh
set -uo pipefail
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
EW=/mnt/data/vcfcache_benchmarks/external_wgs
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
BP37=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh37_joint_af001/blueprint/vcfcache.bcf
VP=/results/vep116/mission/parts
PR=/results/cohort40/prep; mkdir -p $PR/sites
ls /results /results/cohort40 $PR $VP > /dev/null 2>&1   # clears a stale "no such directory" on some nodes
NX=${NX:-1262899}   # records per part of extra (as the new parts of the VEP 116.3 mission)
W=$(mktemp -d /mnt/data/tmp/c40prep-XXXX); trap 'rm -rf $W' EXIT; cd $W
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname)"
test -e $BP && test -e $BP37 && test -e $VP/new.bcf || { echo "inputs missing"; exit 1; }
awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
  $c["status"] == "PASS" && $c["relatedness_status"] == "PASS" {
    print $c["cohort"], $c["sample"], $c["role"], ($c["superpopulation"] == "" ? "NA" : $c["superpopulation"]), $c["assembly"], $c["path"], $c["records"]}' OFS='\t' \
  $EW/manifests/external_wgs_qc.tsv > samples.tsv
awk -F'\t' '{print $1, $3, $5}' samples.tsv | sort | uniq -c
[ $(awk -F'\t' '$5 == "GRCh38" && $3 == "evaluation"' samples.tsv | wc -l) = 40 ] || { echo "not 40 GRCh38 evaluation genomes"; exit 1; }
[ $(awk -F'\t' '$5 == "GRCh38" && $3 == "training"' samples.tsv | wc -l) = 6 ] || { echo "not 6 GRCh38 training genomes"; exit 1; }
cp samples.tsv $PR/samples.tsv
sites() { bcftools view -G $1 -Ou | bcftools norm -m- -Ou 2>/dev/null | bcftools annotate -x INFO,ID,QUAL,FILTER -Ob -W -o $2; }
n=0; while IFS=$'\t' read c s r p a f k; do
  sites $f st_$s.bcf & n=$((n+1)); [ $((n % 8)) = 0 ] && wait
done < samples.tsv; wait
for f in st_*.bcf; do test -s $f.csi || { echo "FAILED sites $f"; exit 1; }; done
awk -F'\t' '$5 == "GRCh38" {print $2}' samples.tsv | while read s; do cp st_$s.bcf $PR/sites/$s.bcf && cp st_$s.bcf.csi $PR/sites/$s.bcf.csi; done
log "sites of $(ls st_*.bcf | wc -l) genomes"
u() { local o=$1; shift; bcftools concat -a -d exact "$@" -Ob -W -o $o --threads 8 2> /dev/null; }
beyond() { bcftools isec -n=1 -w1 -c none $1 $2 -Ob -o $3 2> /dev/null && bcftools index -f $3; }   # records of $1 not in $2
u train.bcf $(awk -F'\t' '$5 == "GRCh38" && $3 == "training" {print "st_" $2 ".bcf"}' samples.tsv)
u all.bcf $(awk -F'\t' '$5 == "GRCh38" {print "st_" $2 ".bcf"}' samples.tsv)
u kpgp.bcf $(awk -F'\t' '$1 == "kpgp" {print "st_" $2 ".bcf"}' samples.tsv)
beyond all.bcf $BP new46.bcf; beyond kpgp.bcf $BP kpgp_new.bcf; beyond new46.bcf $VP/new.bcf extra.bcf
nk=$(bcftools index -n kpgp_new.bcf); nv=$(bcftools index -n $VP/new.bcf); nn=$(bcftools index -n new46.bcf); nx=$(bcftools index -n extra.bcf)
log "46 genomes: $(bcftools index -n all.bcf) sites, $nn beyond the blueprint; KPGP beyond it $nk (VEP mission set $nv); beyond those $nx"
[ $nk = $nv ] && [ $((nv + nx)) = $nn ] || { echo "FAILED: the KPGP set differs from the VEP mission set"; exit 1; }
for f in train new46 extra; do cp $f.bcf $PR/ && cp $f.bcf.csi $PR/; done
# parts by POS (bcftools view -R ... --regions-overlap pos): every record in exactly one part
split_regions() {  # bcf records-per-part prefix
  bcftools query -f '%CHROM\t%POS\n' $1 | awk -v per=$2 -v pre=$3 'BEGIN{OFS="\t"; k=0; n=0}
    { if (n >= per && !($1 == c && $2 == last)) { print c, s, last > (pre "-" k ".regions"); k++; n=0 }   # never cut within one POS
      if (n == 0 || $1 != c) { if (n > 0) print c, s, last > (pre "-" k ".regions"); c=$1; s=$2 }
      last=$2; n++ }
    END { print c, s, last > (pre "-" k ".regions") }'
}
split_regions $PR/extra.bcf $NX $PR/extra
cp $VP/bp-*.regions $VP/new-*.regions $PR/
for s in bp new extra; do echo "$s parts: $(ls $PR/$s-*.regions | wc -l)"; done
# miss rates: sites of each genome beyond the blueprint, and beyond blueprint plus training set
u train37.bcf $(awk -F'\t' '$5 == "GRCh37" && $3 == "training" {print "st_" $2 ".bcf"}' samples.tsv)
printf "cohort\tsample\trole\tsuperpopulation\tassembly\tsites\tbeyond_blueprint\tbeyond_cohort\n" > missrate.tsv
while IFS=$'\t' read c s r p a f k; do
  if [ $a = GRCh38 ]; then b=$BP; t=train.bcf; else b=$BP37; t=train37.bcf; fi
  beyond st_$s.bcf $b m1.bcf; beyond m1.bcf $t m2.bcf
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $c $s $r $p $a $(bcftools index -n st_$s.bcf) $(bcftools index -n m1.bcf) $(bcftools index -n m2.bcf) >> missrate.tsv
done < samples.tsv
column -t missrate.tsv; cp missrate.tsv $PR/
awk -F'\t' 'NR > 1 && $7 == $6 {print "FAILED: no site of " $2 " in its blueprint"; e=1} END {exit e}' missrate.tsv || exit 1
touch $PR/READY; echo DONE
