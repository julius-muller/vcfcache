#!/bin/bash
# The mission pipeline on chr1 of one genome: fastVEP 0.5.0 with eleven known
# sources and two mock custom tracks, uncached against vcfcache with
#   - the gnomAD AF >= 1% and AF >= 0.1% blueprints,
#   - the misses annotated from the same sources stored in narrow blocks
#     (the uncached baseline always uses the sources as fastVEP advises),
#   - a cohort cache (blueprint plus the other KPGP samples), and
#   - a rerun cache that also holds the sample itself (no misses).
# Usage: mission_chr1.sh <deployment> [sample]
set -uo pipefail
DEP=$1; SAMPLE=${2:-KPGP-00128}
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
A=/results/fastvep05/acmg; S=$A/sources
KP=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp
BPD=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints
OLD=/results/fastvep05/caches/fastvep050_clinical_grch38_af001
OUT=/results/fastvep05/cost_model/mission_chr1; mkdir -p $OUT
R=/mnt/data/tmp/mission-chr1; mkdir -p $R/tmp; cd $R
M16=$R/sa16; M12=$R/sa12; mkdir -p $M16 $M12
log() { echo "[$(date -Is)] $*"; }
b() { local s=$(date +%s); $FV sa-build --no-progress --assembly GRCh38 "$@" > build.$$.log 2>&1 || { cat build.$$.log; return 1; }; log "sa-build $1 $2 ${*: -2} $(( $(date +%s) - s ))s"; }
echo "node $(hostname), fastvep $($FV --version)"

# ---- sources ---------------------------------------------------------------
if [ ! -e $R/SOURCES ]; then
  # as fastVEP advises: the existing chr1 stack (dense sources at 65 kb blocks)
  for f in $A/sa/dense/*; do [ "$(basename $f)" = DONE ] || cp -n $f $M16/; done
  # two mock custom tracks: every 93rd base of chr1, three alleles, two scores
  $PY - $FA 1 <<'PYE' | bgzip -@4 -l1 > mock.vcf.gz
import sys
fa, c = sys.argv[1:]
for line in open(fa + ".fai"):
    f = line.split("\t")
    if f[0] == c:
        n, off, bases, width = map(int, f[1:5]); break
out = sys.stdout
out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=%d>\n" % n)
out.write('##INFO=<ID=P,Number=1,Type=Float,Description="mock prediction">\n##INFO=<ID=Q,Number=1,Type=Float,Description="mock confidence">\n')
out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
with open(fa, "rb") as s:
    for pos in range(1, n + 1, 93):
        s.seek(off + (pos - 1) // bases * width + (pos - 1) % bases)
        r = s.read(1).decode().upper()
        if r not in "ACGT":
            continue
        for k, a in enumerate("ACGT"):
            if a != r:
                h = (pos * 48271 + k * 40503) % 1000003
                out.write("chr1\t%d\t.\t%s\t%s\t.\t.\tP=%.3f;Q=%.2f\n" % (pos, r, a, (h % 991) / 991, (h % 97) / 97))
PYE
  log "mock track rows $(zcat mock.vcf.gz | grep -vc '^#')"
  b --source custom_vcf --name TrackA -i mock.vcf.gz -o $M16/tracka --chunk-bits 16 &
  b --source custom_vcf --name TrackB -i mock.vcf.gz -o $M16/trackb --chunk-bits 16 &
  # the same sources in narrow blocks, for the misses of a cached run
  cp -n $A/sa/dense/*.oga $A/sa/dense/*.osi $M12/
  G=$S/gnomad.exomes.v4.1.sites.chr1.vcf.bgz
  b --source gnomad -i $G -o $M12/gnomad_chr1 --chunk-bits 12 &
  ( python3 $A/scripts/extract_gnomad_scores.py --chrom 1 --url $G --out-dir $R/tmp > extract.log 2>&1
    b --source spliceai -i $R/tmp/spliceai_chr1.vcf -o $M12/spliceai_chr1 --chunk-bits 12
    b --source phylop -i $R/tmp/phylop_chr1.tsv -o $M12/phylop_chr1 --chunk-bits 12 ) &
  ( bcftools view -r NC_000001.11 $S/GCF_000001405.40.gz -Oz -o $R/tmp/dbsnp_chr1.vcf.gz --threads 2
    b --source dbsnp -i $R/tmp/dbsnp_chr1.vcf.gz -o $M12/dbsnp_chr1 --chunk-bits 12 ) &
  wait
  b --source custom_vcf --name TrackA -i mock.vcf.gz -o $M12/tracka --chunk-bits 12 &
  b --source custom_vcf --name TrackB -i mock.vcf.gz -o $M12/trackb --chunk-bits 12 &
  b --source clinvar -i $S/clinvar.vcf.gz -o $M12/clinvar --chunk-bits 12 &
  ( python3 -c 'import shutil,sys,zipfile; shutil.copyfileobj(zipfile.ZipFile(sys.argv[1]).open("revel_with_transcript_ids"), sys.stdout.buffer, 1 << 24)' \
      $S/revel-v1.3_all_chromosomes.zip | awk -F, 'NR==1 || $1=="1"' > $R/tmp/revel_chr1.csv
    b --source revel -i $R/tmp/revel_chr1.csv -o $M12/revel_chr1 --chunk-bits 12 ) &
  ( zcat /results/veprs05/data/AlphaMissense_hg38.tsv.gz | awk -F'\t' '/^#/ || $1=="chr1"' | gzip > $R/tmp/am_chr1.tsv.gz
    b --source alphamissense -i $R/tmp/am_chr1.tsv.gz -o $M12/alphamissense_chr1 --chunk-bits 12 ) &
  wait
  rm -rf $R/tmp/*
  n16=$(ls $M16 | wc -l); n12=$(ls $M12 | wc -l)
  ls -la $M16 $M12 | tee $OUT/sources.txt
  [ $n16 = $n12 ] && touch $R/SOURCES || { echo "source sets differ: $n16 vs $n12"; exit 1; }
fi

# ---- inputs and caches -----------------------------------------------------
anno() {  # in.bcf out.bcf sa-dir threads
  bcftools view $1 -Ov | env RAYON_NUM_THREADS=$4 $FV annotate --input - --output - --output-format vcf \
    --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $3 2> anno.err | bcftools view -Ob -o $2 -W --threads 8
}
mkcache() {  # name annotated.bcf sa-dir-for-misses
  local ROOT=$R/caches/$1 C=$R/caches/$1/cache/$1
  rm -rf $ROOT; mkdir -p $C $ROOT/blueprint; cp -r $OLD/workflow $ROOT/
  sed -e "s#--sa-dir /results/richness/sa/clinical_GRCh38#--sa-dir $3#" $OLD/cache/$(basename $OLD)/annotation.yaml > $C/annotation.yaml
  grep -q -- "--sa-dir $3" $C/annotation.yaml || { echo "recipe not rewritten"; exit 1; }
  sed "s#^temp_dir:.*#temp_dir: $R/tmp#" $OLD/cache/$(basename $OLD)/params.snapshot.yaml > $C/params.snapshot.yaml
  ln -f $2 $C/vcfcache_annotated.bcf; ln -f $2.csi $C/vcfcache_annotated.bcf.csi
}
if [ ! -e $R/CACHES ]; then
  bcftools view $KP/$SAMPLE.GRCh38.small_variants.vcf.gz -r chr1 -Ob1 -W --threads 8 -o in.bcf
  for t in af001 af0001; do
    bcftools view -r chr1 $BPD/gnomad_v4.1_GRCh38_joint_$t/blueprint/vcfcache.bcf -Ob -W -o bp_$t.bcf --threads 8
    s=$(date +%s); anno bp_$t.bcf c_$t.bcf $M16 16; log "cache $t: $(bcftools index -n c_$t.bcf) records, $(( $(date +%s) - s ))s, $(du -h c_$t.bcf | cut -f1)"
  done
  # cohort: biallelic sites of the other samples, and of this sample, not in the 1% blueprint
  sites() { bcftools view -r chr1 -G $1 -Ou | bcftools norm -m- -Ou 2>/dev/null | bcftools annotate -x INFO,ID,QUAL,FILTER -Ob -W -o $2; }
  n=0; for f in $KP/*.small_variants.vcf.gz; do
    case $f in *$SAMPLE*) sites $f self.bcf;; *) n=$((n+1)); sites $f o$n.bcf;; esac
  done
  bcftools concat -a -d exact o*.bcf -Ob -W -o others.bcf --threads 8 2>/dev/null
  bcftools isec -n=1 -w1 -c none others.bcf bp_af001.bcf -Ob -o new_cohort.bcf 2>/dev/null; bcftools index -f new_cohort.bcf
  bcftools concat -a -d exact others.bcf self.bcf -Ob -W -o withself.bcf --threads 8 2>/dev/null
  bcftools isec -n=1 -w1 -c none withself.bcf bp_af001.bcf -Ob -o new_rerun.bcf 2>/dev/null; bcftools index -f new_rerun.bcf
  log "cohort: $n other samples, $(bcftools index -n new_cohort.bcf) sites beyond the blueprint; with the sample itself $(bcftools index -n new_rerun.bcf)"
  for k in cohort rerun; do
    anno new_$k.bcf a_$k.bcf $M16 16
    bcftools concat -a c_af001.bcf a_$k.bcf -Ob -W -o c_$k.bcf --threads 8 2>/dev/null
  done
  for k in af001 af0001 cohort rerun; do mkcache $k $R/c_$k.bcf $M16; mkcache ${k}_narrow $R/c_$k.bcf $M12; done
  rm -f o*.bcf o*.bcf.csi; touch $R/CACHES
fi

# ---- runs ------------------------------------------------------------------
N=$(bcftools index -n in.bcf)
digest() { bcftools view -H $1 | awk -F'\t' 'BEGIN{OFS="\t"} {n=split($8,a,";"); asort(a); s=a[1]; for(i=2;i<=n;i++) s=s";"a[i]; $8=s; print}' | md5sum | cut -c1-12; }
vc() {  # arm cache [--uncached]
  local C=$R/caches/$2/cache/$2
  rm -rf st out.bcf*
  local s=$(date +%s.%N)
  env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $C -i in.bcf -o out.bcf \
    -y $C/params.snapshot.yaml --statistics light --stats-dir st -v ${3:-} > log.txt 2>&1 || { echo "FAILED $1"; tail -5 log.txt; return; }
  local e=$(echo "$(date +%s.%N) - $s" | bc)
  local miss=$(grep -ho "missing_variants: [0-9]*" -r st | head -1 | cut -d' ' -f2)
  [ -e dig.$1 ] || digest out.bcf > dig.$1
  printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$N" "${miss:-NA}" "$e" "$(cat dig.$1)" | tee -a $OUT/runs.tsv
}
: > $OUT/runs.tsv
for rep in 1 2 3; do
  vc uncached af001 --uncached
  for k in af001 af001_narrow af0001 af0001_narrow cohort cohort_narrow rerun; do vc $k $k; done
done
# fastVEP alone, for reference
for i in 1 2; do s=$(date +%s.%N); anno in.bcf d.bcf $M16 8; echo "fastvep_alone $(echo "$(date +%s.%N) - $s" | bc)" | tee -a $OUT/runs.tsv; done
echo DONE
