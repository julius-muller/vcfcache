"""Cost-model table and figure for the whole-genome mission benchmarks.

Reads the per-run tables in results/ (3 KPGP samples x 3 repeats each) and writes
  results/cost_model_table.tsv / .md   one row per annotator x stack x cache x miss-step width
  results/cost_model_figure.png / .pdf (needs matplotlib)

Per configuration, with D = direct (uncached) time, R = rerun time (no misses: vcfcache's own
cost of reading the cache and merging), m = miss fraction, C = cached time:
  observed speed-up  D / C
  ideal speed-up     D / (R + m * D)   if every miss cost the annotator what one variant costs it
                                       in direct annotation
  miss cost ratio    (C - R) / (m * D) what one miss actually costs, relative to that
Each value is the median over repeats per sample, then the median over samples.

Usage: python figure.py [--no-figure]
"""

from __future__ import annotations

import statistics as st
import sys
from pathlib import Path

RES = Path(__file__).resolve().parent / "results"
SAMPLES = ["KPGP-00006", "KPGP-00086", "KPGP-00128"]
# (annotator, sources, file prefix, {arm: (cache, miss-step width)})
SETS = [
    ("fastVEP 0.5.0", 11, "missionfn10", None),
    ("fastVEP 0.5.0", 15, "mission2fn10", None),
    ("vep-rs 0.5.0", 11, "veprs_mission", None),
    ("vep-rs 0.5.0", 15, "veprs_mission2", None),
]
ARMS = {
    "af001": ("gnomAD AF>=0.1%", "65 kb"),
    "af001_narrow": ("gnomAD AF>=0.1%", "4 kb"),
    "af001_n10": ("gnomAD AF>=0.1%", "1 kb"),
    "cohort": ("cohort", "65 kb"),
    "cohort_narrow": ("cohort", "4 kb"),
    "cohort_n10": ("cohort", "1 kb"),
    "rerun": ("rerun (all hits)", "-"),
}


def read(path: Path) -> tuple[dict[str, list[float]], dict[str, int], int, set[str]]:
    times: dict[str, list[float]] = {}
    misses: dict[str, int] = {}
    n = 0
    digests = set()
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) != 5:
            continue  # e.g. "fastvep_alone <s>" rows
        arm, nvar, miss, sec, dig = f
        times.setdefault(arm, []).append(float(sec))
        n = int(nvar)
        if miss != "NA":
            misses[arm] = int(miss)
        digests.add(dig)
    return times, misses, n, digests


def rows() -> list[dict]:
    out = []
    for annot, nsrc, prefix, _ in SETS:
        per = {}  # arm -> list of per-sample dicts
        for s in SAMPLES:
            times, misses, n, digests = read(RES / f"{prefix}_wgs_{s}.tsv")
            assert len(digests) == 1, f"{prefix} {s}: outputs differ"
            d = st.median(times["uncached"])
            r = st.median(times["rerun"])
            for arm in ARMS:
                if arm not in times:
                    continue
                c = st.median(times[arm])
                m = misses[arm] / n
                per.setdefault(arm, []).append(
                    dict(D=d, R=r, C=c, m=m, obs=d / c, ideal=d / (r + m * d),
                         ratio=(c - r) / (m * d) if m else float("nan"))
                )
        for arm, lst in per.items():
            cache, width = ARMS[arm]
            if annot.startswith("vep-rs"):
                width = "-"  # vep-rs looks sources up per variant (tabix)
            row = dict(annotator=annot, sources=nsrc, cache=cache, miss_width=width, samples=len(lst))
            for k in ["D", "R", "C", "m", "obs", "ideal", "ratio"]:
                row[k] = st.median(x[k] for x in lst)
            row["obs_min"] = min(x["obs"] for x in lst)
            row["obs_max"] = max(x["obs"] for x in lst)
            out.append(row)
    return out


def write_tables(rs: list[dict]) -> None:
    head = ["annotator", "sources", "cache", "miss_width", "samples", "direct_s", "rerun_s",
            "cached_s", "miss_pct", "speedup", "speedup_min", "speedup_max", "ideal_speedup",
            "share_of_ideal", "miss_cost_ratio"]
    lines = ["\t".join(head)]
    md = ["| Annotator | Sources | Cache | Miss-step blocks | Direct (s) | Cached (s) | Misses | "
          "Speed-up (range) | Ideal | Share of ideal | Cost per miss |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rs:
        ratio = "" if r["ratio"] != r["ratio"] else f"{r['ratio']:.2f}"
        lines.append("\t".join(map(str, [
            r["annotator"], r["sources"], r["cache"], r["miss_width"], r["samples"],
            f"{r['D']:.0f}", f"{r['R']:.1f}", f"{r['C']:.1f}", f"{100 * r['m']:.2f}",
            f"{r['obs']:.2f}", f"{r['obs_min']:.2f}", f"{r['obs_max']:.2f}", f"{r['ideal']:.2f}",
            f"{r['obs'] / r['ideal']:.2f}", ratio])))
        md.append(
            f"| {r['annotator']} | {r['sources']} | {r['cache']} | {r['miss_width']} | {r['D']:.0f} | "
            f"{r['C']:.0f} | {100 * r['m']:.1f}% | {r['obs']:.2f}x ({r['obs_min']:.2f}-{r['obs_max']:.2f}) | "
            f"{r['ideal']:.1f}x | {100 * r['obs'] / r['ideal']:.0f}% | {ratio + 'x' if ratio else '-'} |")
    (RES / "cost_model_table.tsv").write_text("\n".join(lines) + "\n")
    note = (
        "\nWhole genomes (KPGP-00006, -00086, -00128; 4.6-4.8 M variants), 3 repeats each, median over "
        "repeats then over samples; outputs identical to direct annotation in every run. Direct = the "
        "annotator alone as advised (fastVEP sources at 65 kb blocks). Ideal = D / (R + m D): the speed-up "
        "if every miss cost what one variant costs in direct annotation, with R the rerun time (vcfcache's "
        "own cost). Cost per miss = (C - R) / (m D), relative to that.\n")
    (RES / "cost_model_table.md").write_text("\n".join(md) + "\n" + note)


BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def figure(rs: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": INK2,
                         "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                         "axes.linewidth": 0.8})
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1, 1.15]},
                               facecolor=SURFACE)

    # Panel A: observed vs ideal speed-up
    def colour(r):
        if r["annotator"].startswith("vep-rs"):
            return AQUA
        # fastVEP rerun: no misses, the cache is built from 65 kb sources
        return {"65 kb": BLUE, "1 kb": ORANGE, "-": BLUE}.get(r["miss_width"], GREY)

    marker = {"gnomAD AF>=0.1%": "o", "cohort": "s", "rerun (all hits)": "^"}
    for r in rs:
        if r["miss_width"] == "4 kb":
            continue
        filled = r["sources"] == 15
        c = colour(r)
        a.scatter(r["ideal"], r["obs"], s=46, marker=marker[r["cache"]], zorder=3,
                  facecolor=c if filled else SURFACE, edgecolor=c, linewidth=1.6)
    lim = (1, 200)
    a.plot(lim, lim, color=INK2, lw=1, ls=(0, (4, 3)), zorder=1)
    a.text(60, 75, "every miss costs\none direct variant", color=INK2, fontsize=7.5, rotation=0,
           ha="right", va="bottom")
    a.set_xscale("log"), a.set_yscale("log")
    a.set_xlim(*lim), a.set_ylim(1, 200)
    ticks = [1, 2, 3, 5, 10, 20, 50, 100, 200]
    a.set_xticks(ticks, [f"{t}x" for t in ticks]), a.set_yticks(ticks, [f"{t}x" for t in ticks])
    a.minorticks_off()
    a.axhline(3, color=GRID, lw=1, zorder=0)
    a.text(1.08, 3.1, "3x", color=INK2, fontsize=7.5, va="bottom")
    a.grid(True, color=GRID, lw=0.6, zorder=0)
    a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    a.set_xlabel("Ideal speed-up  D / (R + m·D)")
    a.set_ylabel("Observed speed-up  D / C")
    a.set_title("A  Observed vs ideal speed-up, whole genomes", loc="left", fontsize=9.5,
                color=INK, fontweight="bold")
    # selective direct labels: the fastVEP 15-source points and vep-rs 15-source points
    for r in rs:
        if r["sources"] != 15 or r["miss_width"] == "4 kb" or r["cache"].startswith("rerun"):
            continue
        gn = r["cache"].startswith("gnomAD")
        lab = f"{'0.1%' if gn else 'cohort'}, {r['obs']:.1f}x"
        # cohort labels right of the point; 0.1% labels below-right (clear of the diagonal),
        # or left where a cohort point sits below-right
        if not gn:
            off, ha, va = (9, 0), "left", "center"
        elif r["obs"] < 2.5:
            off, ha, va = (-9, 0), "right", "center"
        else:
            off, ha, va = (6, -7), "left", "top"
        a.annotate(lab, (r["ideal"], r["obs"]), xytext=off, textcoords="offset points",
                   fontsize=7, color=INK2, ha=ha, va=va)
    h = [Line2D([], [], ls="", marker="o", ms=7, mfc=BLUE, mec=BLUE, label="fastVEP, sources at 65 kb"),
         Line2D([], [], ls="", marker="o", ms=7, mfc=ORANGE, mec=ORANGE, label="fastVEP, misses at 1 kb"),
         Line2D([], [], ls="", marker="o", ms=7, mfc=AQUA, mec=AQUA, label="vep-rs"),
         Line2D([], [], ls="", marker="o", ms=7, mfc=INK2, mec=INK2, label="15 sources"),
         Line2D([], [], ls="", marker="o", ms=7, mfc=SURFACE, mec=INK2, mew=1.6, label="11 sources"),
         Line2D([], [], ls="", marker="o", ms=7, mfc=SURFACE, mec=INK2, label="gnomAD AF≥0.1% cache"),
         Line2D([], [], ls="", marker="s", ms=7, mfc=SURFACE, mec=INK2, label="cohort cache"),
         Line2D([], [], ls="", marker="^", ms=7, mfc=SURFACE, mec=INK2, label="rerun (all hits)")]
    a.legend(handles=h, loc="upper left", frameon=False, fontsize=7.2, handletextpad=0.3,
             borderaxespad=0.2, labelcolor=INK)

    # Panel B: where the time goes, fastVEP with 15 sources
    sel = [r for r in rs if r["annotator"].startswith("fastVEP") and r["sources"] == 15]
    d = sel[0]["D"]
    bars = [("Direct (fastVEP alone)", None, d, BLUE, None)]
    for cache, cname in [("gnomAD AF>=0.1%", "gnomAD AF≥0.1% cache"), ("cohort", "cohort cache")]:
        for w in ["65 kb", "4 kb", "1 kb"]:
            r = next(x for x in sel if x["cache"] == cache and x["miss_width"] == w)
            bars.append((f"{cname}, misses at {w}", r["R"], r["C"] - r["R"],
                         {"65 kb": BLUE, "1 kb": ORANGE}.get(w, GREY), r["obs"]))
    r = next(x for x in sel if x["cache"].startswith("rerun"))
    bars.append(("Rerun (all hits)", r["R"], 0, None, r["obs"]))
    y = list(range(len(bars)))[::-1]
    for yi, (lab, rr, ann, c, obs) in zip(y, bars):
        left = 0
        if rr:
            b.barh(yi, rr, color="#c9c8c2", height=0.62, edgecolor=SURFACE, linewidth=2)
            left = rr
        if ann:
            b.barh(yi, ann, left=left, color=c, height=0.62, edgecolor=SURFACE, linewidth=2)
        total = (rr or 0) + ann
        sp = "" if obs is None else f"   {obs:.2f}x"
        b.text(total + 12, yi, f"{total:.0f} s{sp}", va="center", fontsize=7.5, color=INK)
    b.set_yticks(y, [x[0] for x in bars], fontsize=7.8, color=INK)
    b.tick_params(axis="y", length=0)
    b.set_xlim(0, d * 1.22)
    b.set_xlabel("Wall time per genome (s), 8 cores")
    b.grid(True, axis="x", color=GRID, lw=0.6)
    b.set_axisbelow(True)
    for s in ("top", "right", "left"):
        b.spines[s].set_visible(False)
    b.set_title("B  fastVEP, 15 sources: where the time goes", loc="left",
                fontsize=9.5, color=INK, fontweight="bold")
    hb = [plt.Rectangle((0, 0), 1, 1, color="#c9c8c2", label="vcfcache's own cost (rerun time)"),
          plt.Rectangle((0, 0), 1, 1, color=BLUE, label="fastVEP, sources at 65 kb"),
          plt.Rectangle((0, 0), 1, 1, color=GREY, label="fastVEP, sources at 4 kb"),
          plt.Rectangle((0, 0), 1, 1, color=ORANGE, label="fastVEP, sources at 1 kb")]
    b.legend(handles=hb, loc="lower right", frameon=False, fontsize=7.2, labelcolor=INK)
    for ax in (a, b):
        ax.set_facecolor(SURFACE)
    fig.tight_layout(w_pad=2.5)
    for ext in ("png", "pdf"):
        fig.savefig(RES / f"cost_model_figure.{ext}", dpi=300, facecolor=SURFACE)


if __name__ == "__main__":
    rs = rows()
    write_tables(rs)
    print((RES / "cost_model_table.md").read_text())
    if "--no-figure" not in sys.argv:
        figure(rs)
