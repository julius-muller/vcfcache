#!/bin/bash
#SBATCH --job-name=integration-confirm
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=01:00:00
# Two narrowly targeted controls, using the original frozen deployments.
# No uncached baseline is rerun: the original direct outputs are audited.
set -euo pipefail
kind=$1
head_root=/results/deployment/vcfcache-eval-head-63757a07
main_root=/results/deployment/vcfcache-eval-main-136bafa9
python=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
output=/results/branch-integration/confirm-${kind}-${SLURM_JOB_ID}
work=/mnt/data/tmp/integration-confirm-${kind}-${SLURM_JOB_ID}
mkdir -p "$output" "$work"
case "$kind" in
  fastvep)
    original=/mnt/data/tmp/branch-eval-fastvep-1300
    cache=$original/cache/fastvep_clinical_grch38_af001
    scale=wes
    # A small timing control, not a new technical-replicate campaign.
    options=(--arms auto native --paired-order --repeats 2) ;;
  vep)
    original=/mnt/data/tmp/branch-eval-vep-1317
    cache=$original/cache/vep115.2_everything
    scale=wgs
    # Reverse the original main -> legacy@head order to test the apparent gap.
    options=(--arms legacy main-legacy --legacy-root "$main_root" --repeats 1) ;;
  *) echo "unknown control: $kind" >&2; exit 2 ;;
esac
test -f "$cache/vcfcache_annotated.bcf"
test -f "$original/runs/${scale}-direct-1/output.bcf"
{ echo "kind=$kind"; echo "node=$(hostname)"; echo "job=$SLURM_JOB_ID";
  echo "source_head=63757a07"; echo "source_main=136bafa9"; } > "$output/CONDITIONS.txt"
cd "$work"
PYTHONPATH="$head_root" "$python" -m benchmarks.adaptive_routing \
  --cache "$cache" --input "/results/branch-eval/inputs/$scale.bcf" \
  --output "$output/benchmark" --work-dir "$work/runs" --cpus 16 \
  --slurm-steps --reference-root "$original/runs" "${options[@]}"
