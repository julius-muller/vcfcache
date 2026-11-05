#!/usr/bin/env bash
# Private, reproducible toolchain; never replaces /usr/local/bin.
set -euo pipefail
prefix=$1
mkdir -p "$prefix"
test ! -e "$prefix/READY"
cd "$prefix"
curl -fL --retry 3 -o bcftools-1.25.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.25/bcftools-1.25.tar.bz2
sha256sum bcftools-1.25.tar.bz2 > SOURCE.sha256
tar -xjf bcftools-1.25.tar.bz2
cd bcftools-1.25
# Match the deployed zlib-based 1.21 codec configuration: no compression A/B.
./configure --prefix="$prefix" --without-libdeflate
make -j8
make install
make -C htslib-1.25 prefix="$prefix" install
"$prefix/bin/bcftools" --version > "$prefix/VERSION.txt"
"$prefix/bin/bcftools" --version-only | grep -Fx '1.25+htslib-1.25'
sha256sum "$prefix/bin/bcftools" > "$prefix/BINARY.sha256"
touch "$prefix/READY"
