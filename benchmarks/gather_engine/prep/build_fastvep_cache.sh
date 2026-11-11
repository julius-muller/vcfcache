#!/bin/bash
# Build the fastVEP $V clinical vcfcache cache for one assembly.
set -euo pipefail
ASM=$1; asm=${ASM,,}; V=${2:?fastVEP version, e.g. 0.5.0}; TAG=fastvep${V//[.-]/}
D=/mnt/data/apps/fastvep/$V
test -s $D/node_ready.json
if [ $ASM = GRCh38 ]; then BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf; FA=$D/data/GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa
else BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh37_joint_af001/blueprint/vcfcache.bcf; FA=$D/data/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa; fi
TC=$D/data/$ASM/Homo_sapiens.$ASM.fastvep.cache
SA=/results/richness/sa/clinical_$ASM
OLD=/results/richness/caches/fastvep_clinical_${asm}_af001/cache/fastvep_clinical_${asm}_af001
NAME=${TAG}_clinical_${asm}_af001
ROOT=/results/fastvep05/caches/$NAME; C=$ROOT/cache/$NAME
test ! -e $ROOT
mkdir -p $C/auxiliary $ROOT/blueprint $ROOT/workflow
cp /results/native-bcf/veprs-caches-20261023b/veprs_core/workflow/init.yaml $ROOT/workflow/init.yaml
cp $OLD/annotation.yaml $C/annotation.yaml
grep -q -- "--sa-dir $SA" $C/annotation.yaml
cat > $C/params.snapshot.yaml <<Y
genome_build: $ASM
annotation_tool_cmd: $D/bin/fastvep
tool_version_command: $D/bin/fastvep --version
bcftools_cmd: bcftools
temp_dir: /mnt/data/tmp
threads: 8
rayon_threads: 8
transcript_cache: $TC
fasta: $FA
optional_checks: {}
Y
start=$(date +%s)
bcftools view $BP -Ov | env RAYON_NUM_THREADS=16 $D/bin/fastvep annotate --input - --output - --output-format vcf \
  --transcript-cache $TC --fasta $FA --hgvs --no-progress --sa-dir $SA 2> $C/auxiliary/fastvep_stderr.txt | \
  bcftools view -Ob -o $C/vcfcache_annotated.bcf -W --threads 8
seconds=$(( $(date +%s) - start ))
in=$(bcftools index -n $BP); out=$(bcftools index -n $C/vcfcache_annotated.bcf); echo "in=$in out=$out seconds=$seconds"
test "$in" = "$out"
for t in CSQ FV_CLINVAR FV_GNOMAD; do bcftools view -h $C/vcfcache_annotated.bcf | grep -q "##INFO=<ID=$t," ; done
cat > $C/build_provenance.json <<J
{"kind": "fastvep_${V}_clinical_from_bundled_blueprint", "assembly": "$ASM", "tool_version": "$($D/bin/fastvep --version)",
 "build_seconds": $seconds, "records": $out, "node": "$(hostname)", "supplementary": "$SA (built with fastVEP 0.3.0)",
 "cache_sha256": "$(sha256sum $C/vcfcache_annotated.bcf | cut -d' ' -f1)"}
J
cat $C/build_provenance.json; ls -la $C
touch /results/fastvep05/caches/$NAME.COMPLETE
