#!/bin/bash
#SBATCH --job-name=branch-eval-prepare
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=01:00:00
# Stage one sample at four nested scales plus two startup probes.
# Every consumer of a bcftools pipe drains it; an early-exiting consumer would
# raise SIGPIPE under `set -o pipefail` and abort this job with a false failure.
set -euo pipefail
sample=LP6005519-DNA_A06
source=/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/sgdp/${sample}.GRCh38.small_variants.vcf.gz
panel_bed=/results/richness/regions/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed
wes_bed=/results/richness/regions/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed
out=/results/branch-eval/inputs
work=/mnt/data/tmp/branch-eval-prepare-${SLURM_JOB_ID}
mkdir -p "$out" "$work"

bcftools view "$source" -Ob -W --threads 8 -o "$out/wgs.bcf"
bcftools view -R "$wes_bed" "$out/wgs.bcf" -Ob -W --threads 8 -o "$out/wes.bcf"
bcftools view -R "$panel_bed" "$out/wgs.bcf" -Ob -W --threads 8 -o "$out/panel.bcf"

# One locus per occupied contig: the smallest input that still forces the
# engine to open every contig a real exome touches.
bcftools view -H "$out/wes.bcf" | awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' > "$work/minimal.tsv"
sort -k1,1 -k2,2n "$work/minimal.tsv" > "$work/minimal.sorted.tsv"
bcftools view -T "$work/minimal.sorted.tsv" "$out/wes.bcf" -Ob -W -o "$out/minimal.bcf"
sed -n '1p' "$work/minimal.sorted.tsv" > "$work/single.tsv"
bcftools view -T "$work/single.tsv" "$out/wes.bcf" -Ob -W -o "$out/single.bcf"

cp "$work/minimal.sorted.tsv" "$work/single.tsv" "$out/"
{
  printf 'input\trecords\tcontigs\n'
  for name in single minimal panel wes wgs; do
    records=$(bcftools index -n "$out/$name.bcf")
    contigs=$(bcftools view -H "$out/$name.bcf" | awk -F'\t' '!seen[$1]++ {n++} END {print n+0}')
    printf '%s\t%s\t%s\n' "$name" "$records" "$contigs"
  done
} > "$out/INVENTORY.tsv"
sha256sum "$out"/*.bcf "$panel_bed" "$wes_bed" > "$out/SOURCE_SHA256SUMS"
cat "$out/INVENTORY.tsv"
