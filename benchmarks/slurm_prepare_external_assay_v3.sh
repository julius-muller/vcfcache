#!/usr/bin/env bash
#SBATCH --partition=benchmark
#SBATCH --cpus-per-task=4
#SBATCH --mem=8000M
#SBATCH --time=01:00:00

set -euo pipefail

manifest=${VCFCACHE_ASSAY_PREP_MANIFEST:?required}
task_id=${SLURM_ARRAY_TASK_ID:?required}
row=$(awk -F '\t' -v id="$task_id" 'NR == id + 2 {print; exit}' "$manifest")
[[ -n "$row" ]]
IFS=$'\t' read -r observed_id cohort sample assay source_vcf bed output_vcf <<<"$row"
[[ "$observed_id" == "$task_id" ]]
[[ -s "$source_vcf" && -s "$source_vcf.tbi" && -s "$bed" ]]

if [[ -s "$output_vcf" && -s "$output_vcf.tbi" ]]; then
  bcftools index --nrecords "$output_vcf"
  exit 0
fi

mkdir -p "$(dirname "$output_vcf")"
partial="${output_vcf%.vcf.gz}.partial-$SLURM_JOB_ID.vcf.gz"
trap 'rm -f -- "$partial" "$partial.tbi"' EXIT
bcftools view --regions-file "$bed" --output-type z --threads 4 \
  --output "$partial" "$source_vcf"
bcftools index --tbi --threads 4 "$partial"
records=$(bcftools index --nrecords "$partial")
((records > 0))
mv "$partial" "$output_vcf"
mv "$partial.tbi" "$output_vcf.tbi"
python3 - "$output_vcf.provenance.json" "$source_vcf" "$bed" "$records" \
  "$cohort" "$sample" "$assay" <<'PY'
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

target, source, bed = map(Path, sys.argv[1:4])
document = {
    "created_at": datetime.now(timezone.utc).isoformat(),
    "source_vcf": str(source),
    "source_vcf_sha256": digest(source),
    "interval_bed": str(bed),
    "interval_bed_sha256": digest(bed),
    "output_records": int(sys.argv[4]),
    "cohort": sys.argv[5],
    "sample": sys.argv[6],
    "assay": sys.argv[7],
    "command": "bcftools view --regions-file BED --output-type z",
}
target.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
PY
trap - EXIT
