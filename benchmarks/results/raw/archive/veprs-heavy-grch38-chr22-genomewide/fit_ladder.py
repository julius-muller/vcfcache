# Two-point fit per rung: t = s + c*N, using the MEDIAN of the three reps at
# each size so a cold first touch cannot set the slope on its own.
import csv, statistics as st

rows = list(csv.DictReader(open("heavy_1251.tsv"), delimiter="\t"))
by = {}
for r in rows:
    if r["rc"] != "0":
        continue
    by.setdefault(r["config"], {}).setdefault(r["size"], []).append(
        (float(r["wall_seconds"]), int(r["records"])))

THR_LEGACY = 27.5e-6 / 0.9258        # 29.7 us/variant
THR_FUSED  = (27.5e-6 / 1.986) / 0.9258   # 15.0 us/variant, from the fused check

print(f"{'rung':26s} {'c (us/var)':>10s} {'s (s)':>7s} {'rec/s':>8s}  {'legacy':>7s} {'fused':>7s}")
print("-" * 72)
fits = {}
for cfg, sizes in by.items():
    if "small" not in sizes or "big" not in sizes:
        continue
    t1 = st.median([t for t, _ in sizes["small"]]); n1 = sizes["small"][0][1]
    t2 = st.median([t for t, _ in sizes["big"]]);   n2 = sizes["big"][0][1]
    c = (t2 - t1) / (n2 - n1)
    s = t1 - n1 * c
    fits[cfg] = c
    v = lambda thr: "WIN" if c > thr else "lose"
    print(f"{cfg:26s} {c*1e6:10.2f} {s:7.3f} {n2/t2:8.0f}  {v(THR_LEGACY):>7s} {v(THR_FUSED):>7s}")

print("\nper-source cost of a whole-file genome-wide gnomAD lookup:")
base = fits["L0_core"]
for k, n in [("L1_gw_af", 1), ("W2_gw_wholefile_x2", 2),
             ("W4_gw_wholefile_x4", 4), ("W8_gw_wholefile_x8", 8)]:
    if k in fits:
        print(f"  {n} instance(s): {fits[k]*1e6:7.2f} us total, "
              f"{(fits[k]-base)/n*1e6:6.2f} us per source")

print("\ngenome-wide vs exome-restricted, one whole-file source:")
for k, label in [("L1_gw_af", "genome-wide joint (99.9 M rec)"),
                 ("X1_exome_wholefile", "exome-restricted (183.7 M rec)")]:
    if k in fits:
        print(f"  {label:34s} {fits[k]*1e6:6.2f} us/variant")
