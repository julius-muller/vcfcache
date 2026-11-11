#!/bin/bash
# Stage fastVEP $V and the Ensembl 116 GFF3 on shared storage.
set -euo pipefail
V=${1:?fastVEP version, e.g. 0.5.0}
ST=/results/apps-staging/fastvep; mkdir -p $ST; cd $ST
T=fastVEP-v$V-x86_64-unknown-linux-gnu.tar.gz
U=https://github.com/Huang-lab/fastVEP/releases/download/v$V
[ -s $T ] || curl -fsSL -o $T $U/$T
curl -fsSL -o $T.sha256 $U/$T.sha256
sha256sum -c $T.sha256
G=Homo_sapiens.GRCh38.116.gff3.gz
[ -s $G ] || { curl -fsSL -o $G.part https://ftp.ensembl.org/pub/release-116/gff3/homo_sapiens/$G && mv $G.part $G; }
gzip -t $G; sha256sum $G > $G.sha256
touch $ST/DOWNLOADED.$V; ls -la $ST
