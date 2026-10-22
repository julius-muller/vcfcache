# Two-point fit of VCFcache cached-workflow cost, legacy vs fused engine.
# t = s + c*N, from the repo's own benchmarks.fused_engine harness.
N1, N2 = 200_000, 800_000
runs = {"legacy": (2.9585575889941538, 9.577666146011325),
        "fused":  (1.839520955007174, 5.171885628995369)}
F = 0.9258            # measured median WGS hit rate
PROD_C = 27.5e-6      # production per-variant overhead, legacy (analyze_breakeven.py)

fit = {}
print(f"{'engine':8s} {'c (us/var)':>11s} {'s (s)':>8s}")
for k, (t1, t2) in runs.items():
    c = (t2 - t1) / (N2 - N1)
    s = t1 - N1 * c
    fit[k] = (c, s)
    print(f"{k:8s} {c*1e6:11.3f} {s:8.3f}")

ratio = fit["legacy"][0] / fit["fused"][0]
print(f"\nper-variant overhead reduction: {ratio:.3f}x")
print(f"start-up barely moves: {fit['legacy'][1]:.3f} s -> {fit['fused'][1]:.3f} s "
      f"({fit['fused'][1]/fit['legacy'][1]:.3f}x) -- consistent with a per-variant-only win")

old_thr = PROD_C / F
new_c   = PROD_C / ratio
new_thr = new_c / F
print(f"\nthreshold = (VCFcache per-variant overhead) / f")
print(f"  legacy: {PROD_C*1e6:.1f} us / {F} = {old_thr*1e6:5.1f} us/variant")
print(f"  fused : {new_c*1e6:.2f} us / {F} = {new_thr*1e6:5.1f} us/variant  (scaled by the measured ratio)")

print(f"\nwhere the vep-rs rungs land:")
for name, c in [("core (job 1242, head)", 5.68e-6),
                ("core (job 1250, worker)", 7.49e-6),
                ("--everything (job 1242)", 14.00e-6),
                ("+genome-wide gnomAD (job 1250 L1)", 34.38e-6)]:
    def verdict(thr):
        return "cache WINS" if c > thr else "cache loses"
    print(f"  {name:36s} {c*1e6:6.2f} us   legacy: {verdict(old_thr):10s}   fused: {verdict(new_thr)}")
