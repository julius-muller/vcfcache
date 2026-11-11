"""Print the BED of part K of N: contiguous contig runs of comparable record count."""

import subprocess
import sys

bp, n, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rows = [
    line.split("\t")
    for line in subprocess.check_output(
        ["bcftools", "index", "--stats", bp], text=True
    ).splitlines()
]
occ = [
    (r[0], int(r[1]) if r[1].isdigit() else 2**31 - 1, int(r[2]))
    for r in rows
    if len(r) >= 3 and int(r[2])
]
target = sum(r[2] for r in occ) / n
parts, cur, load = [], [], 0
for i, (c, length, rec) in enumerate(occ):
    cur.append((c, length))
    load += rec
    if load >= target and i < len(occ) - 1 and len(parts) + 1 < n:
        parts.append(cur)
        cur, load = [], 0
parts.append(cur)
assert len(parts) == n, (len(parts), n)
for c, length in parts[k]:
    print(f"{c}\t0\t{length}")
