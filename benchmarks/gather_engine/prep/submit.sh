#!/bin/bash
# Submit the whole preparation chain. Run once on the controller as appuser from /tmp.
set -euo pipefail
P=/results/apps-staging/prep; L=/results/apps-staging/logs; mkdir -p /mnt/data/slurm-results/apps-staging/logs
sb() { sbatch --parsable --partition=benchmark --chdir=/tmp --output=$L/%x-%j.out "$@" | cut -d';' -f1; }
NODES=(sl-w1 sl-w2 sl-w3 sl-w4 sl-w5 sl-w6); N=${#NODES[@]}
d38=$(sb --job-name=dl-vep-GRCh38 --nodelist=sl-w1 -c2 --mem=4000M --time=12:00:00 --wrap "bash $P/dl_vep.sh GRCh38")
d37=$(sb --job-name=dl-vep-GRCh37 --nodelist=sl-w2 -c2 --mem=4000M --time=12:00:00 --wrap "bash $P/dl_vep.sh GRCh37")
dfv=$(sb --job-name=dl-fastvep --nodelist=sl-w3 -c1 --mem=2000M --time=02:00:00 --wrap "bash $P/dl_fastvep.sh")
echo "downloads: $d38 $d37 $dfv"
p38=(); p37=(); fv=()
for k in $(seq 0 $((N-1))); do n=${NODES[$k]}
  f=$(sb --job-name=inst-fastvep-$n --nodelist=$n -c8 --mem=16000M --time=02:00:00 --dependency=afterok:$dfv --wrap "bash $P/install_fastvep.sh")
  i38=$(sb --job-name=inst-vep38-$n --nodelist=$n -c4 --mem=12000M --time=06:00:00 --dependency=afterok:$d38 --wrap "bash $P/install_vep.sh GRCh38")
  b38=$(sb --job-name=vep38-part$k --nodelist=$n --exclusive -c16 --mem=30000M --time=10:00:00 --dependency=afterok:$i38 --wrap "bash $P/build_vep_part.sh GRCh38 $k $N")
  i37=$(sb --job-name=inst-vep37-$n --nodelist=$n -c4 --mem=12000M --time=06:00:00 --dependency=afterok:$d37:$i38 --wrap "bash $P/install_vep.sh GRCh37")
  b37=$(sb --job-name=vep37-part$k --nodelist=$n --exclusive -c16 --mem=30000M --time=10:00:00 --dependency=afterok:$i37,afterany:$b38 --wrap "bash $P/build_vep_part.sh GRCh37 $k $N")
  p38+=($b38); p37+=($b37); fv+=($f)
  echo "$n fastvep=$f inst38=$i38 part38=$b38 inst37=$i37 part37=$b37"
done
a38=$(sb --job-name=vep38-assemble --nodelist=sl-w1 -c16 --mem=20000M --time=03:00:00 --dependency=afterok:$(IFS=:; echo "${p38[*]}") --wrap "bash $P/assemble_vep.sh GRCh38 $N")
a37=$(sb --job-name=vep37-assemble --nodelist=sl-w2 -c16 --mem=20000M --time=03:00:00 --dependency=afterok:$(IFS=:; echo "${p37[*]}") --wrap "bash $P/assemble_vep.sh GRCh37 $N")
c38=$(sb --job-name=fastvep-cache-GRCh38 --nodelist=sl-w3 -c16 --mem=28000M --time=06:00:00 --dependency=afterok:${fv[2]} --wrap "bash $P/build_fastvep_cache.sh GRCh38")
c37=$(sb --job-name=fastvep-cache-GRCh37 --nodelist=sl-w4 -c16 --mem=28000M --time=06:00:00 --dependency=afterok:${fv[3]} --wrap "bash $P/build_fastvep_cache.sh GRCh37")
echo "assemble: $a38 $a37  fastvep caches: $c38 $c37"
