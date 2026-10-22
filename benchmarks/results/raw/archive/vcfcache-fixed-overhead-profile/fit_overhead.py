# Medians per cell, then the fixed/per-variant split.
# Arm C of job 1253 is EXCLUDED: its three inputs were all single-contig
# regardless of the span they claimed (the genome-wide pool is contig-ordered,
# so "first 250 matching rows" always landed in chr1). Job 1254 supersedes it.
import csv, statistics as st

def load(path, keyf):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["rc"] != "0" or not r["wall_seconds"]:
            continue
        out.setdefault(keyf(r), []).append(float(r["wall_seconds"]))
    return {k: st.median(v) for k, v in out.items()}

A = load("overhead_1253/overhead_1253.tsv",
         lambda r: (r["arm"], r["cell"], r["engine"], int(r["records"])))
S = load("span_1254/span_1254.tsv",
         lambda r: (r["engine"], int(r["contigs"]), int(r["records"])))

print("ARM A -- size at genome-wide span (22 contigs; gw1 is 1 contig), zero misses")
print(f"{'N':>7} {'legacy':>9} {'fused':>9}")
sizes = [1, 250, 2000, 75000]
for n in sizes:
    l = A.get(("A", f"gw{n}", "legacy", n)); f = A.get(("A", f"gw{n}", "fused", n))
    print(f"{n:>7} {l:>9.1f} {f:>9.1f}")

l250, l75k = A[("A","gw250","legacy",250)], A[("A","gw75000","legacy",75000)]
f250, f75k = A[("A","gw250","fused",250)],  A[("A","gw75000","fused",75000)]
for name, a, b in (("legacy", l250, l75k), ("fused", f250, f75k)):
    c = (b - a) / (75000 - 250); s = a - 250 * c
    print(f"\n{name}: fixed {s:6.1f} s, per-variant {c*1e6:5.1f} us")
    for n in (250, 75000):
        tot = s + n*c
        print(f"   N={n:>6}: fixed is {100*s/tot:4.1f}% of VCFcache's cost")

print("\n\nSPAN SWEEP (job 1254) -- 250 variants held constant, span varied")
print(f"{'contigs':>8} {'legacy':>9} {'fused':>9}")
for k in (1, 2, 4, 8, 22):
    rows = [(e, v) for (e, kk, _), v in S.items() if kk == k for e in [e]]
    l = next(v for (e, kk, _), v in S.items() if kk == k and e == "legacy")
    f = next(v for (e, kk, _), v in S.items() if kk == k and e == "fused")
    print(f"{k:>8} {l:>9.1f} {f:>9.1f}")

print("\n\nARM B -- real assays, with misses, so VEP runs")
print(f"{'assay':<10} {'N':>7} {'miss':>5} {'legacy':>8} {'fused':>8} {'VCFcache':>9} {'VEP':>7}")
for cell, n, miss, zero in [("panel", 203, 3, l250), ("wes", 75135, 1339, l75k)]:
    pre = [k for k in A if k[0] == "B" and k[1].startswith(cell) and k[2] == "legacy" and k[3] == n]
    pf  = [k for k in A if k[0] == "B" and k[1].startswith(cell) and k[2] == "fused"  and k[3] == n]
    lg, fu = A[pre[0]], A[pf[0]]
    print(f"{cell:<10} {n:>7} {miss:>5} {lg:>8.1f} {fu:>8.1f} {zero:>9.1f} {lg-zero:>7.1f}")
