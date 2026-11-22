#!/bin/bash
# Resubmit the VEP 116.3 chain after the downloads (install, build, assemble, remove 115, repoint vep-rs).
set -euo pipefail
P=/results/apps-staging/prep; L=/results/apps-staging/logs
sb() { sbatch --parsable --partition=benchmark --chdir=/tmp --output=$L/%x-%j.out "$@" | cut -d';' -f1; }
NODES=(sl-w1 sl-w2 sl-w3 sl-w4 sl-w5 sl-w6); N=${#NODES[@]}
p38=(); p37=(); rm=()
for k in $(seq 0 $((N-1))); do n=${NODES[$k]}
  i38=$(sb --job-name=inst-vep38-$n --nodelist=$n -c4 --mem=12000M --time=06:00:00 --wrap "bash $P/install_vep.sh GRCh38")
  b38=$(sb --job-name=vep38-part$k --nodelist=$n --exclusive -c16 --mem=30000M --time=10:00:00 --dependency=afterok:$i38 --wrap "bash $P/build_vep_part.sh GRCh38 $k $N")
  i37=$(sb --job-name=inst-vep37-$n --nodelist=$n -c4 --mem=12000M --time=06:00:00 --dependency=afterok:$i38 --wrap "bash $P/install_vep.sh GRCh37")
  b37=$(sb --job-name=vep37-part$k --nodelist=$n --exclusive -c16 --mem=30000M --time=10:00:00 --dependency=afterok:$i37,afterany:$b38 --wrap "bash $P/build_vep_part.sh GRCh37 $k $N")
  r=$(sb --job-name=rm-vep115-$n --nodelist=$n -c1 --mem=2000M --time=01:00:00 --dependency=afterok:$i37 --wrap "bash $P/remove_vep115.sh")
  p38+=($b38); p37+=($b37); rm+=($r)
  echo "$n inst38=$i38 part38=$b38 inst37=$i37 part37=$b37 rm115=$r"
done
a38=$(sb --job-name=vep38-assemble --nodelist=sl-w1 -c16 --mem=20000M --time=03:00:00 --dependency=afterok:$(IFS=:; echo "${p38[*]}") --wrap "bash $P/assemble_vep.sh GRCh38 $N")
a37=$(sb --job-name=vep37-assemble --nodelist=sl-w2 -c16 --mem=20000M --time=03:00:00 --dependency=afterok:$(IFS=:; echo "${p37[*]}") --wrap "bash $P/assemble_vep.sh GRCh37 $N")
rp=$(sb --job-name=repoint-veprs -c1 --mem=1000M --time=00:10:00 --dependency=afterok:$(IFS=:; echo "${rm[*]}") --wrap "bash $P/repoint_veprs.sh")
echo "assemble: $a38 $a37 repoint: $rp"
