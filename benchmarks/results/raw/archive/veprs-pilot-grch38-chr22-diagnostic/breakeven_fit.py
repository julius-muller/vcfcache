rows = [
 ("core", 0.346339430, 2.038518157),
 ("symbol_canonical", 0.344428845, 2.241170297),
 ("hgvs", 0.374169915, 3.528315094),
 ("everything", 0.393701621, 4.565737449),
]
N1, N2 = 16385, 314314
dN = N2 - N1
VC_S, VC_C = 21.7, 27.5e-6   # vcfcache overhead: s + c*N (fastVEP fit)
F = 0.9258                    # median WGS hit rate

print(f"{'config':18s} {'c (us/var)':>11s} {'s (s)':>7s} {'rec/s':>9s} {'f*c-27.5':>10s} {'needed c':>9s}")
for name, t1, t2 in rows:
    c = (t2 - t1) / dN
    s = t1 - N1 * c
    print(f"{name:18s} {c*1e6:11.2f} {s:7.3f} {N2/t2:9.0f} {(F*c-VC_C)*1e6:10.2f} {VC_C/F*1e6:9.1f}")

print()
print("WGS projection, N = 4.8e6, everything config (richest reachable):")
c = (4.565737449 - 0.393701621) / dN
s = 0.393701621 - N1 * c
N = 4.8e6
td = s + N*c
tc = s + (1-F)*N*c + VC_S + VC_C*N
print(f"  direct  {td:8.1f} s   ({td/60:.1f} min)")
print(f"  cached  {tc:8.1f} s   ({tc/60:.1f} min)")
print(f"  speedup {td/tc:8.2f}x")
