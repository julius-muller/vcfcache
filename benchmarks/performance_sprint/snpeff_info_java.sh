#!/usr/bin/env bash
# Clinical recipe contract: annotate INFO, preserve the sample's original ID.
# Match the historical cache's ASCII text-encoding contract explicitly;
# Java 21 otherwise defaults to UTF-8 regardless of LC_ALL=C.
# Apply the supported SnpSift switch to both database stages, not to SnpEff.
set -euo pipefail
args=()
sift=false
for argument in "$@"; do
    if [[ "$argument" == */SnpSift.jar ]]; then
        sift=true
    fi
    args+=("$argument")
    if $sift && [[ "$argument" == annotate ]]; then
        args+=(-noId)
    fi
done
exec /results/snpeff/dist/jre/bin/java -Dfile.encoding=US-ASCII "${args[@]}"
