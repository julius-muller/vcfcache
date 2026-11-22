#!/bin/bash
# Join the annotated parts into a vcfcache cache directory for VEP 116.3.
set -euo pipefail
ASM=$1; N=$2
if [ "$ASM" = GRCh38 ]; then BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf
else BP=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh37_joint_af001/blueprint/vcfcache.bcf; fi
OLD=/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_${ASM}_joint_af001/cache/vep115.2_everything
IN=/results/vep116/build-$ASM
ROOT=/results/vep116/caches/gnomad_v4.1_${ASM}_joint_af001
C=$ROOT/cache/vep116.3_everything
test ! -e $C
mkdir -p $C $ROOT/blueprint $ROOT/workflow
cp /results/native-bcf/veprs-caches-20261023b/veprs_core/workflow/init.yaml $ROOT/workflow/init.yaml
parts=$(for k in $(seq 0 $((N-1))); do echo $IN/part-$k.bcf; done)
bcftools concat $parts -Ob --threads 16 -o $C/vcfcache_annotated.bcf
bcftools index -f --threads 16 $C/vcfcache_annotated.bcf
in=$(bcftools index -n $BP); out=$(bcftools index -n $C/vcfcache_annotated.bcf); ref=$(bcftools index -n $OLD/vcfcache_annotated.bcf)
echo "blueprint=$in cache=$out vep115_cache=$ref"
# VEP skips records on contigs it has no annotation for, as 115 did: the keys must equal the VEP 115 cache.
test "$out" = "$ref"
cmp <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $OLD/vcfcache_annotated.bcf | md5sum) <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $C/vcfcache_annotated.bcf | md5sum)
echo "variant keys identical to the VEP 115 cache"
sed -e 's/115\.2/116.3/g' $OLD/annotation.yaml > $C/annotation.yaml
grep -q 'required_tool_version: "116.3"' $C/annotation.yaml
cat > $C/params.snapshot.yaml <<Y
genome_build: "$ASM"
annotation_tool_cmd: "apptainer exec -B /mnt/data:/mnt/data /mnt/data/apps/ensembl-vep/116/vep.sif vep"
tool_version_command: "apptainer exec -B /mnt/data:/mnt/data /mnt/data/apps/ensembl-vep/116/vep.sif vep | awk '/ensembl-vep/ {print \$NF}'"
threads: 8
bcftools_cmd: "bcftools"
temp_dir: "/tmp"
vep_cache: "/mnt/data/apps/ensembl-vep/116/cachedir"
vep_buffer: 100000
vep_forks: 8
Y
python3 - <<P > $C/build_provenance.json
import json, glob, subprocess
parts = [json.load(open(f)) for f in sorted(glob.glob("$IN/part-*.json"))]
print(json.dumps(dict(kind="vep116.3_everything_from_bundled_blueprint", assembly="$ASM", parts=parts,
    node_seconds=sum(p["seconds"] for p in parts), records=$out, blueprint_records=$in, vep115_cache_records=$ref, forks_per_part=10,
    vep=open("/results/apps-staging/ensembl-vep-116/vep.version").read().strip(),
    cache_sha256=subprocess.check_output(["sha256sum", "$C/vcfcache_annotated.bcf"], text=True).split()[0]), indent=1))
P
cat $C/build_provenance.json; ls -la $C
bcftools view -h $C/vcfcache_annotated.bcf | grep -c "##INFO=<ID=CSQ,"
touch /results/vep116/caches/$ASM.COMPLETE
