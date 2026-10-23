#!/bin/bash
#SBATCH --job-name=native-bcf-transcripts
#SBATCH --partition=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=02:00:00
set -euo pipefail
destination=$1
test ! -e "$destination"
mkdir -p "$destination"
veprs=/results/veprs
source=/mnt/data/apps/ensembl-vep/115/cachedir/homo_sapiens/115_GRCh38
apptainer exec -B "$source:/input:ro" -B "$destination:/output" \
  -B "$veprs/dist:/scripts:ro" "$veprs/dist/ensembl-vep_116.2.sif" \
  perl /scripts/storable_to_json.pl /input /output 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y \
  > "$destination/conversion.out" 2> "$destination/conversion.err"
test -s "$destination/info.json"
for chromosome in {1..22} X Y; do
  test -d "$destination/transcripts/$chromosome"
done
# Converter reports failures as warnings; do not mistake an incomplete cache for
# a successful conversion even when its process exit status is zero.
test ! -s "$destination/conversion.err"
sha256sum "$source/info.txt" "$veprs/dist/storable_to_json.pl" "$veprs/dist/vep" \
  > "$destination/SOURCE_SHA256SUMS"
touch "$destination/COMPLETE"
