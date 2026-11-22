#!/bin/bash
# Smoke run, gate, then the full cohort, all from one deployment.
# Usage: launch_cohort.sh <deployment> <commit> <results root>
set -euo pipefail
DEP=$1; COMMIT=$2; ROOT=$3
PY=/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python
L=/results/apps-staging/logs
run() { env PYTHONPATH=$DEP $PY -m benchmarks.gather_engine.cohort "$@"; }
run plan $ROOT-smoke --deployment $DEP --commit $COMMIT --smoke
smoke=$(python3 -c "import json;print(json.load(open('$ROOT-smoke/submissions.json'))['jobs']['smoke'])")
gate=$(sbatch --parsable --partition=benchmark --chdir=/tmp -c1 --mem=2000M --time=00:20:00 \
  --dependency=afterok:$smoke --job-name=cohort-gate --output=$L/%x-%j.out \
  --wrap "PYTHONPATH=$DEP $PY -m benchmarks.gather_engine.cohort gate $ROOT-smoke" | cut -d';' -f1)
echo "smoke=$smoke gate=$gate"
run plan $ROOT --deployment $DEP --commit $COMMIT --after $gate
