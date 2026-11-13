#!/bin/bash
# Dense-source probe: does fastVEP's cost for a source that scores every
# possible SNV (SpliceAI2 / CADD shape) grow with the chunks touched or with
# the variants annotated? One chromosome is the miniature: variants and cache
# misses per chunk are the same as genome-wide, so times scale with length.
#
# Mock sources (deterministic, no licence): "spliceai" through fastVEP's own
# SpliceAI parser, and "mockdb", ten numeric INFO columns through custom_vcf.
# Usage: run.sh <engine deployment> [START END] (default: all of chr21)
set -euo pipefail
DEP=$1; START=${2:-1}; END=${3:-46709983}; TAG=${4:-chr21}
C=21; REG=chr$C:$START-$END
export PATH=/results/input-preserving-20261105-tools/bin:$PATH
FVD=/mnt/data/apps/fastvep/0.5.0; FV=$FVD/bin/fastvep
TC=$FVD/data/GRCh38/Homo_sapiens.GRCh38.fastvep.cache; FA=$FVD/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
SRC=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp/KPGP-00128.GRCh38.small_variants.vcf.gz
BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
OLD=/results/fastvep05/acmg/caches/fastvep050_B0_chr1
R=/mnt/data/tmp/dense-probe-$TAG; mkdir -p $R/sa $R/caches; cd $R   # node-local disk
OUT=/results/fastvep05/dense/$TAG; mkdir -p $OUT
log() { echo "[$(date -Is)] $*"; }
echo "node $(hostname), $(nproc) cpus, region $REG, fastvep $($FV --version), disk $(df --output=source $R | tail -1)"

# --- mock sources -----------------------------------------------------------
seq_lines() {  # the chromosome's sequence lines, 60 bases each
  $PY - $FA $C <<'PY'
import sys, shutil
fa, c = sys.argv[1:]
for line in open(fa + ".fai"):
    f = line.split("\t")
    if f[0] == c:
        n, off, bases, width = map(int, f[1:5]); break
with open(fa, "rb") as s:
    s.seek(off); left = n + n // bases * (width - bases) + width
    while left > 0:
        b = s.read(min(left, 1 << 24)); left -= len(b); sys.stdout.buffer.write(b)
        if not b: break
PY
}
mock() {  # kind: spliceai | mockdb
  seq_lines | awk -v kind=$1 -v s=$START -v e=$END 'BEGIN{split("A C G T",B," ");
      print "##fileformat=VCFv4.2";
      if(kind=="spliceai") print "##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"mock\">";
      else for(j=1;j<=10;j++) print "##INFO=<ID=S" j ",Number=1,Type=Float,Description=\"mock\">";
      print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}
    { n=length($0); for(i=1;i<=n;i++){ pos++; if(pos<s) continue; if(pos>e) exit;
        r=toupper(substr($0,i,1)); if(r=="N") continue;
        for(k=1;k<=4;k++){ a=B[k]; if(a==r) continue; h=(pos*48271+k*40503)%1000003;
          if(kind=="spliceai")
            # nine in ten delta scores are 0.00, as in gene bodies away from splice sites
            printf "chr21\t%d\t.\t%s\t%s\t.\t.\tSpliceAI=%s|G%d|%.2f|%.2f|%.2f|%.2f|%d|%d|%d|%d\n", pos,r,a,a,int(pos/40000),
              (h%10==0)*(h%97)/100,(h%11==0)*(h%89)/100,(h%13==0)*(h%83)/100,(h%17==0)*(h%79)/100,
              h%101-50,h%97-48,h%89-44,h%83-41;
          else
            printf "chr21\t%d\t.\t%s\t%s\t.\t.\tS1=%.3f;S2=%.3f;S3=%.3f;S4=%.3f;S5=%.3f;S6=%.2f;S7=%.2f;S8=%.2f;S9=%d;S10=%d\n", pos,r,a,
              (h%991)/991,(h%983)/983,(h%977)/977,(h%971)/971,(h%967)/967,(h%953)/9.53,(h%947)/9.47,(h%941)/94.1,h%37,h%1009;
        } } }' | bgzip -@6 -l1 > $R/$1.vcf.gz
}
if [ ! -e $R/MOCKS ]; then
  mock spliceai & mock mockdb & wait
  for k in spliceai mockdb; do
    n=$(zcat $R/$k.vcf.gz | grep -vc '^#'); echo "{\"source\":\"$k\",\"records\":$n,\"vcf_gz_bytes\":$(stat -c %s $R/$k.vcf.gz)}" | tee $OUT/mock_$k.json
  done
  touch $R/MOCKS
fi
mkdir -p $R/sa/base; ln -f /results/fastvep05/acmg/sa/stack/clinvar.osa2 $R/sa/base/ 2>/dev/null || cp /results/fastvep05/acmg/sa/stack/clinvar.osa2 $R/sa/base/
# Each source is built once per chunk width; the stacks link to the result.
sab() {  # source name bits
  local o=$R/sa/build/$2_$3; [ -e $o.osa2 ] && return
  local t=$(date +%s); /usr/bin/time -f "sa-build $2 bits $3: %e s, peak %M KB" $FV sa-build --no-progress --assembly GRCh38 --source $1 \
    $([ $1 = custom_vcf ] && echo --name $2) -i $R/$2.vcf.gz -o $R/sa/build/tmp_$2_$3 --chunk-bits $3 && mv $R/sa/build/tmp_$2_$3.osa2 $o.osa2
}
mkdir -p $R/sa/build
# custom_vcf builds in memory (about 170 bytes a record: 120 M records were
# killed at 20 GB), so the region is sized for it; the SpliceAI parser streams.
for bits in 16 12; do sab spliceai spliceai $bits & sab custom_vcf mockdb $bits & wait; test -e $R/sa/build/mockdb_$bits.osa2 -a -e $R/sa/build/spliceai_$bits.osa2; done
for bits in 16 12; do for st in D1 D2; do
  S=$R/sa/${st}_$bits; mkdir -p $S; cp -l -f $R/sa/base/clinvar.osa2 $S/; cp -l -f $R/sa/build/spliceai_$bits.osa2 $S/spliceai.osa2
  [ $st = D1 ] || cp -l -f $R/sa/build/mockdb_$bits.osa2 $S/mockdb.osa2
done; done
ls -la $R/sa/D*/ $R/sa/base | tee $OUT/sa_files.txt

# --- inputs -----------------------------------------------------------------
bcftools view -r $REG $SRC -Ob1 -W --threads 8 -o full.bcf
bcftools view -r $REG $BP -Ob -W --threads 8 -o bp.bcf
log "input $(bcftools index -n full.bcf) records, blueprint $(bcftools index -n bp.bcf)"
cat > p.yaml <<Y
genome_build: GRCh38
annotation_tool_cmd: $FV
tool_version_command: $FV --version
bcftools_cmd: bcftools
temp_dir: $R
threads: 8
rayon_threads: 8
transcript_cache: $TC
fasta: $FA
optional_checks: {}
Y
fv() { env RAYON_NUM_THREADS=8 $FV annotate --transcript-cache $TC --fasta $FA --hgvs --no-progress --output-format vcf "$@"; }
build_cache() {  # stack
  local NAME=dense_$1 ROOT=$R/caches/dense_$1; local K=$ROOT/cache/$NAME
  [ -e $ROOT/COMPLETE ] && return
  rm -rf $ROOT; mkdir -p $K $ROOT/blueprint; cp -r $OLD/workflow $ROOT/
  sed -e "s#--sa-dir /results/fastvep05/acmg/sa/stack#--sa-dir $R/sa/$1#" $OLD/cache/$(basename $OLD)/annotation.yaml > $K/annotation.yaml
  grep -q -- "--sa-dir $R/sa/$1" $K/annotation.yaml; cp p.yaml $K/params.snapshot.yaml
  bcftools view bp.bcf -Ov | fv --input - --output - --sa-dir $R/sa/$1 2> $ROOT/err.txt | bcftools view -Ob -o $K/vcfcache_annotated.bcf -W --threads 8
  test "$(bcftools index -n bp.bcf)" = "$(bcftools index -n $K/vcfcache_annotated.bcf)"; touch $ROOT/COMPLETE
}
STACKS="base D1_16 D1_12 D2_16 D2_12"
for s in $STACKS; do build_cache $s; done
vc() {  # stack out [--uncached]
  rm -rf st $2 $2.csi; env PYTHONPATH=$DEP $PY -c "from vcfcache.cli import main; main()" annotate -a $R/caches/dense_$1/cache/dense_$1 \
    -i full.bcf -o $2 -y p.yaml --statistics light --stats-dir st -v --debug ${3:-} > log_vc.txt 2>&1
}
# The misses of a cached run are the same records for every stack.
vc base warm.bcf; cp $(find st -name misses.bcf | head -1) misses.bcf; bcftools index -f misses.bcf
bcftools view full.bcf -Ov -o full.vcf; bcftools view misses.bcf -Ov -o misses.vcf
log "misses $(bcftools index -n misses.bcf)"
$PY - > $OUT/chunks.json <<PY
import json, subprocess
out = {}
for name in ("full", "misses"):
    pos = [int(l.split()[1]) for l in open(f"$R/{name}.vcf") if l[0] != "#"]
    out[name] = {"records": len(pos), **{f"chunks_{b}": len({p >> b for p in pos}) for b in (20, 16, 12)}}
out["region"] = {f"chunks_{b}": ($END >> b) - ($START >> b) + 1 for b in (20, 16, 12)}
print(json.dumps(out, indent=1))
PY
cat $OUT/chunks.json

# --- A: fastVEP alone, full input against misses only -------------------------
evict() { $PY -c 'import os,sys
for p in sys.argv[1:]:
    fd=os.open(p,os.O_RDONLY); os.posix_fadvise(fd,0,0,os.POSIX_FADV_DONTNEED); os.close(fd)' "$@"; }
: > $OUT/direct.tsv
for rep in 1 2; do for s in $STACKS; do for inp in full misses; do for state in cold warm; do
  [ $state = cold ] && { [ $rep = 1 ] || continue; evict $R/sa/$s/*.osa2; }
  /usr/bin/time -v -o tm.txt env RAYON_NUM_THREADS=8 $FV annotate --transcript-cache $TC --fasta $FA --hgvs --no-progress \
     --output-format vcf --input $inp.vcf --output o.vcf --sa-dir $R/sa/$s 2> err.txt
  g() { grep "$1" tm.txt | awk '{print $NF}'; }
  echo -e "$s\t$inp\t$state\t$rep\t$(g 'Elapsed (wall')\t$(g 'User time')\t$(g 'System time')\t$(g 'Maximum resident')\t$(g 'File system inputs')\t$(grep -vc '^#' o.vcf)\t$(grep -c 'SpliceAI=' o.vcf)\t$(grep -c 'FV_MOCKDB=' o.vcf || true)\t$(grep -v '^#' o.vcf | awk 'length($4)==1 && length($5)==1 && $5!="*"' | wc -l)" | tee -a $OUT/direct.tsv
done; done; done; done
grep -m3 "FV_MOCKDB=" o.vcf | cut -c1-300 > $OUT/example_records.txt; grep "^##INFO" o.vcf | grep -iv "ID=CSQ" | cut -c1-200 >> $OUT/example_records.txt

# --- B: vcfcache end to end, identical output required -------------------------
: > $OUT/vcfcache.tsv
for rep in 1 2 3; do for s in $STACKS; do for mode in cached --uncached; do
  t=$(date +%s.%N); vc $s out.bcf $([ $mode = cached ] || echo $mode); e=$(echo "$(date +%s.%N) - $t" | bc)
  # Second digest: INFO tags sorted within a record and records sorted, since
  # bcftools annotate can swap a hit and a miss that share a position.
  echo -e "$s\t$mode\t$rep\t$e\t$(bcftools view -H out.bcf | md5sum | cut -c1-12)\t$(bcftools view -H out.bcf | awk -F'\t' -v OFS='\t' '{n=split($8,t,";"); asort(t); x=t[1]; for(i=2;i<=n;i++) x=x ";" t[i]; $8=x; print}' | LC_ALL=C sort | md5sum | cut -c1-12)" | tee -a $OUT/vcfcache.tsv
done; done; done
log DONE
