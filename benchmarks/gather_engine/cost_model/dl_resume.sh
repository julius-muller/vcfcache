#!/bin/bash
# Resume the dbNSFP 5.4a download after a broken transfer and check its md5. Usage: dl_resume.sh
set -uo pipefail
ls /results > /dev/null 2>&1
D=/results/sources/dbnsfp; U=$(cat $D/URL_GRCH38)
for try in 1 2 3 4 5 6; do
  curl -fsSL --retry 20 --retry-delay 30 -C - -o $D/dbNSFP5.4a_grch38.gz $U && break; echo "try $try ended, $(du -h $D/dbNSFP5.4a_grch38.gz | cut -f1)"; sleep 30
done
a=$(md5sum $D/dbNSFP5.4a_grch38.gz | cut -d' ' -f1); grep -q $a $D/dbNSFP5.4a_grch38.gz.md5 || { echo "md5 mismatch"; exit 1; }
echo "dbNSFP md5 ok, $(du -h $D/dbNSFP5.4a_grch38.gz | cut -f1)"; touch /results/sources/DOWNLOADED
