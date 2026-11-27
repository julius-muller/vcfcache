"""Speed-up per assay and annotator, and the cost of building a cache: tables and two figures.

Reads the per-run tables in results/ (3 KPGP samples x 3 repeats each; exome = the genome restricted
to the Twist core exome pad 125, panel = ACMG SF v3.3 CDS pad 20) and writes
  results/assay_table.tsv / .md         annotator x assay x cache: direct, cached, speed-up
  results/build_cost_table.tsv / .md    annotator: cache build time and break-even sample count
  results/assay_figure.png / .pdf       speed-up per assay (needs matplotlib)
  results/build_cost_figure.png / .pdf  time per sample including the cache build, and the turnaround
                                        of one sample with a published cache (needs matplotlib)
Speed-up = direct / cached, median over repeats per sample, then over samples (range = samples).
Every annotator runs its mission stack of 11 sources (SnpEff: the 8 VCF sources of the Ensembl VEP
stack). fastVEP's miss step uses 1 kb sources for genomes and exomes and its 65 kb sources for panels
(the rule from assay_test.sh). Outputs are checked: one digest per assay over all arms (fastVEP,
vep-rs), or record-by-record comparisons that show only the known annotator quirks (Ensembl VEP:
HGNC_ID, CLIN_SIG order; SnpEff: INFO tags in header order, snpeff_info_order.sh).

The build cost B is the time to annotate the 18.9 M-record gnomAD AF>=1% blueprint (one node, 16
threads; results/cache_build_times.tsv, SnpEff from results/snpeff_mission_parts.tsv). A user who
builds the cache for their own stack and annotates N samples spends B/N + C per sample instead of D;
break-even N* = B / (D - C). With a published cache B = 0.

Usage: python assay_figure.py [--no-figure]
"""

from __future__ import annotations

import math
import statistics as st
import sys
from pathlib import Path

RES = Path(__file__).resolve().parent / "results"
SAMPLES = ["KPGP-00006", "KPGP-00086", "KPGP-00128"]
ASSAYS = [("wgs", "Genome"), ("wes", "Exome"), ("panel", "Panel")]
CACHES = [("af001", "gnomAD AF≥1%"), ("cohort", "Cohort"), ("rerun", "Rerun")]
N10 = {"af001": "af001_n10", "cohort": "cohort_n10", "rerun": "rerun"}
ANNOTATORS = [
    # name, genome runs, assay runs, genome compare, assay compare, allowed comparison lines, arm map
    ("fastVEP 0.5.0", "missionfn10_wgs_{}.tsv", "assay_{}.tsv", None, None, None,
     {"wgs": N10, "wes": N10}),
    ("vep-rs 0.5.0", "veprs_mission_wgs_{}.tsv", "veprs_assay_{}.tsv", None, None, None, {}),
    ("Ensembl VEP 116.3", "vep116_mission_wgs_{}.tsv", "vep116_assay_{}.tsv", "vep116_compare_wgs_{}.tsv",
     "vep116_assay_compare_{}.tsv",
     ("records", "identical", "HGNC_ID only", "HGNC_ID genes", "HGNC_ID and list order only",
      "list order only"), {}),
    ("SnpEff 5.4c", "snpeff_mission_wgs_{}.tsv", "snpeff_assay_{}.tsv", "snpeff_compare_wgs_{}.tsv",
     "snpeff_assay_compare_{}.tsv",
     ("records", "identical", "identical (INFO tag order only)", "ANN warnings only", "ANN warnings"), {}),
]


def runs(path: Path, assay: str) -> list[tuple[str, int, str, float, str]]:
    """(arm, variants, misses, seconds, digest) rows of one assay from a runs table."""
    out = []
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if line.startswith("#") or f[0] == "assay":
            continue
        if f[0] in ("wes", "panel"):  # fastVEP assay_test.sh: assay, arm, variants, misses, s, digest
            a, f = f[0], f[1:]
        elif "." in f[0]:  # <assay>.<arm>
            a, arm = f[0].split(".", 1)
            f = [arm] + f[1:]
        else:
            a = "wgs"
        if a != assay or len(f) < 5:
            continue
        out.append((f[0], int(f[1]), f[2], float(f[3]), f[4]))
    return out


def check_compare(path: Path, allowed: tuple[str, ...], assay: str, arms: set[str]) -> None:
    """Assert that a record-by-record comparison shows only the allowed categories for every arm."""
    seen = set()
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if assay != "wgs":
            if f[0] != assay:
                continue
            f = f[1:]
        if len(f) != 3 or f[0] not in arms:  # skips the timing rows some comparison files start with
            continue
        assert f[1].startswith(allowed), f"{path.name} {assay}: {f[0]} {f[1]}"
        seen.add(f[0])
    assert seen == arms, f"{path.name} {assay}: compared {sorted(seen)}, need {sorted(arms)}"


def rows() -> list[dict]:
    """One row per annotator x assay x cache (medians over samples)."""
    out = []
    for name, gfile, afile, gcmp, acmp, allowed, armmap in ANNOTATORS:
        for assay, _ in ASSAYS:
            per: dict[str, list[dict]] = {}
            for s in SAMPLES:
                path = RES / (gfile if assay == "wgs" else afile).format(s)
                if not path.exists():
                    break
                rs = [r for r in runs(path, assay) if r[0] != "warmup"]
                times: dict[str, list[float]] = {}
                miss, n = {}, 0
                for arm, nv, m, sec, _ in rs:
                    times.setdefault(arm, []).append(sec)
                    n = nv
                    if m != "NA":
                        miss[arm] = int(m)
                amap = armmap.get(assay, {k: k for k, _ in CACHES})
                if allowed is None:
                    digs = {r[4] for r in runs(path, assay) if r[0] in set(amap.values()) | {"uncached", "warmup"}}
                    assert len(digs) == 1, f"{name} {assay} {s}: outputs differ {digs}"
                else:
                    cmp_path = RES / (gcmp if assay == "wgs" else acmp).format(s)
                    check_compare(cmp_path, allowed, assay, set(amap.values()))
                d = st.median(times["uncached"])
                for cache, _ in CACHES:
                    c = st.median(times[amap[cache]])
                    per.setdefault(cache, []).append(
                        dict(D=d, C=c, n=n, m=miss[amap[cache]] / n, sp=d / c, reps=len(times[amap[cache]])))
            for cache, label in CACHES:
                lst = per.get(cache, [])
                if len(lst) != len(SAMPLES):
                    continue
                out.append(dict(annotator=name, assay=assay, cache=label, samples=len(lst),
                                variants=st.median(x["n"] for x in lst), D=st.median(x["D"] for x in lst),
                                C=st.median(x["C"] for x in lst), m=st.median(x["m"] for x in lst),
                                sp=st.median(x["sp"] for x in lst), sp_min=min(x["sp"] for x in lst),
                                sp_max=max(x["sp"] for x in lst)))
    return out


def build_costs() -> dict[str, float]:
    """Seconds to annotate the AF>=1% blueprint per annotator (11-source stacks; median of builds)."""
    b: dict[str, list[float]] = {}
    for line in (RES / "cache_build_times.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        if f[1] == "11":
            b.setdefault(f[0], []).append(float(f[3]))
    parts = RES / "snpeff_mission_parts.tsv"
    if parts.exists():  # part, records, seconds, node
        b["SnpEff 5.4c"] = [sum(float(x.split("\t")[2]) for x in parts.read_text().splitlines()
                                if x.startswith("bp-"))]
    return {k: st.median(v) for k, v in b.items()}


def write_tables(rs: list[dict], b: dict[str, float]) -> list[dict]:
    """Write the assay and build-cost tables; return the build-cost rows."""
    names = dict(ASSAYS)
    head = ["annotator", "assay", "cache", "samples", "variants", "direct_s", "cached_s", "miss_pct",
            "speedup", "speedup_min", "speedup_max"]
    tsv = ["\t".join(head)]
    md = ["| Annotator | Assay | Cache | Variants | Direct (s) | Cached (s) | Misses | Speed-up (range) |",
          "|---|---|---|---|---|---|---|---|"]
    for r in rs:
        tsv.append("\t".join(map(str, [r["annotator"], r["assay"], r["cache"], r["samples"], int(r["variants"]),
                                       f"{r['D']:.1f}", f"{r['C']:.1f}", f"{100 * r['m']:.2f}", f"{r['sp']:.2f}",
                                       f"{r['sp_min']:.2f}", f"{r['sp_max']:.2f}"])))
        md.append(f"| {r['annotator']} | {names[r['assay']]} | {r['cache']} | {int(r['variants']):,} | "
                  f"{r['D']:.0f} | {r['C']:.1f} | {100 * r['m']:.1f}% | "
                  f"{r['sp']:.2f}x ({r['sp_min']:.2f}-{r['sp_max']:.2f}) |")
    (RES / "assay_table.tsv").write_text("\n".join(tsv) + "\n")
    (RES / "assay_table.md").write_text("\n".join(md) + "\n\nMedian over 3 repeats per sample, then over "
                                        "3 KPGP samples (range over samples). Direct = the annotator alone.\n")
    brows = []
    for name, *_ in ANNOTATORS:
        if name not in b:
            continue
        row = dict(annotator=name, B=b[name])
        for assay in ("wgs", "wes"):
            r = next((x for x in rs if x["annotator"] == name and x["assay"] == assay
                      and x["cache"].startswith("gnomAD")), None)
            if r:
                row[assay] = (r["D"], r["C"], b[name] / (r["D"] - r["C"]))
        brows.append(row)
    tsv = ["annotator\tbuild_s\tgenome_direct_s\tgenome_cached_s\tgenome_break_even\t"
           "exome_direct_s\texome_cached_s\texome_break_even"]
    md = ["| Annotator | Cache build (h) | Genome: direct / cached (s) | Break-even (genomes) | "
          "Exome: direct / cached (s) | Break-even (exomes) |", "|---|---|---|---|---|---|"]
    for r in brows:
        g, e = r.get("wgs"), r.get("wes")
        cells = [f"{x[0]:.0f}\t{x[1]:.0f}\t{x[2]:.1f}" if x else "\t\t" for x in (g, e)]
        tsv.append(f"{r['annotator']}\t{r['B']:.0f}\t{cells[0]}\t{cells[1]}")
        mcell = [f"{x[0]:.0f} / {x[1]:.0f} | {x[2]:.1f}" if x else "- | -" for x in (g, e)]
        md.append(f"| {r['annotator']} | {r['B'] / 3600:.2f} | {mcell[0]} | {mcell[1]} |")
    (RES / "build_cost_table.tsv").write_text("\n".join(tsv) + "\n")
    (RES / "build_cost_table.md").write_text(
        "\n".join(md) + "\n\nCache build = annotating the 18.9 M-record gnomAD AF>=1% blueprint with the "
        "mission stack on one node (16 threads; Ensembl VEP and SnpEff: the sum of 12 parts, one node "
        "each). "
        "Break-even = build / (direct - cached): the number of samples after which building the cache "
        "for one's own stack has paid for itself. With a published cache there is no build. fastVEP's "
        "1 kb miss-step sources are a one-time build of their own, not included.\n")
    return brows


SHORT = {"fastVEP 0.5.0": "fastVEP", "vep-rs 0.5.0": "vep-rs", "Ensembl VEP 116.3": "Ensembl VEP",
         "SnpEff 5.4c": "SnpEff"}
BLUE, ORANGE, AQUA, VIOLET ="#257de0", "#f47230", "#00b894", "#6a4ce2"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def style(plt) -> None:
    """Shared matplotlib settings (as figure.py)."""
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": INK2,
                         "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                         "axes.linewidth": 0.8})


def tidy(ax) -> None:
    """Recessive grid and axes."""
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis="y", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def assay_figure(rs: list[dict]) -> None:
    """Speed-up per assay: one panel per assay, annotators on x, one dot per cache type."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    style(plt)
    names = [a[0] for a in ANNOTATORS if any(r["annotator"] == a[0] for r in rs)]
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), sharey=True, facecolor=SURFACE)
    colours = {"af001": BLUE, "cohort": ORANGE, "rerun": AQUA}
    ticks = [1, 2, 3, 5, 10, 20, 50, 100, 200]
    for ax, (assay, title) in zip(axes, ASSAYS, strict=True):
        for i, name in enumerate(names):
            for j, (cache, label) in enumerate(CACHES):
                r = next((x for x in rs if x["annotator"] == name and x["assay"] == assay
                          and x["cache"] == label), None)
                if r is None:
                    continue
                x = i + (j - 1) * 0.22
                ax.plot([x, x], [r["sp_min"], r["sp_max"]], color=colours[cache], lw=2, zorder=2,
                        solid_capstyle="round")
                ax.scatter(x, r["sp"], s=48, color=colours[cache], edgecolor=SURFACE, linewidth=1.5, zorder=3)
                if cache == "af001":
                    ax.annotate(f"{r['sp']:.1f}x", (x, r["sp"]), xytext=(-6, 0), textcoords="offset points",
                                ha="right", va="center", fontsize=7, color=INK2)
        ax.axhline(1, color=INK2, lw=0.9, zorder=1)
        ax.set_yscale("log")
        ax.set_ylim(0.8, 400)
        ax.set_yticks(ticks, [f"{t}x" for t in ticks])
        ax.minorticks_off()
        labels = []
        for name in names:
            d = next((x["D"] for x in rs if x["annotator"] == name and x["assay"] == assay), None)
            labels.append(SHORT[name] + ("" if d is None else f"\ndirect {minutes(d)}"))
        ax.set_xticks(range(len(names)), labels, fontsize=7.8, color=INK)
        ax.tick_params(axis="x", length=0)
        ax.set_xlim(-0.6, len(names) - 0.4)
        tidy(ax)
        v = next((x["variants"] for x in rs if x["assay"] == assay), 0)
        ax.set_title(f"{title}  ({v / 1e6:.1f} M variants)" if v > 1e6 else f"{title}  ({v:,.0f} variants)",
                     loc="left", fontsize=9.5, color=INK, fontweight="bold")
    axes[0].set_ylabel("Speed-up over direct annotation (log)")
    h = [Line2D([], [], ls="", marker="o", ms=7, color=colours[c], label=lab) for c, lab in
         [("af001", "gnomAD AF≥1% cache"), ("cohort", "cohort cache (+22 KPGP samples)"),
          ("rerun", "rerun (sample already cached)")]]
    axes[2].legend(handles=h, loc="upper right", frameon=False, fontsize=7.2, labelcolor=INK)
    fig.text(0.01, 0.01, "Dot: median of 3 samples (3 repeats each); bar: range over samples.",
             fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1), w_pad=1.5)
    for ext in ("png", "pdf"):
        fig.savefig(RES / f"assay_figure.{ext}", dpi=300, facecolor=SURFACE)


def build_figure(brows: list[dict], rs: list[dict]) -> None:
    """Time per sample with the cache build over the number of samples (A, B); turnaround (C).

    C: one sample with a published cache (no build), direct against the two cached arms.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    style(plt)
    colours = dict(zip([a[0] for a in ANNOTATORS], [BLUE, ORANGE, AQUA, VIOLET], strict=True))
    fig = plt.figure(figsize=(15, 4.3), facecolor=SURFACE, layout="constrained")
    fig.get_layout_engine().set(rect=(0, 0.04, 1, 0.96), wspace=0.04)
    grid = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.1])
    a0 = fig.add_subplot(grid[0])
    axes = [a0, fig.add_subplot(grid[1], sharey=a0)]
    axes[1].tick_params(labelleft=False)
    ns = [1.2 ** k for k in range(0, 41)]  # 1 .. ~1500 samples
    for ax, (assay, title, nmax) in zip(axes, [("wgs", "Genomes", 50), ("wes", "Exomes", 1000)], strict=True):
        xs = [n for n in ns if n <= nmax]
        ends = []
        for r in brows:
            if assay not in r:
                continue
            d, c, be = r[assay]
            col = colours[r["annotator"]]
            ys = [(r["B"] / n + c) / d for n in xs]
            ax.plot(xs, ys, color=col, lw=2, zorder=3)
            if be <= nmax:
                ax.scatter(be, 1, s=48, color=col, edgecolor=SURFACE, linewidth=1.5, zorder=4)
            ends.append([math.log10(ys[-1]), f"{SHORT[r['annotator']]}  (break-even {be:.1f})"])
        ends.sort()
        for k in range(1, len(ends)):  # keep end labels apart (log10 units)
            ends[k][0] = max(ends[k][0], ends[k - 1][0] + 0.09)
        for y, lab in ends:
            ax.text(xs[-1] * 1.08, 10 ** y, lab, fontsize=7, color=INK2, va="center")
        ax.axhline(1, color=INK2, lw=0.9, zorder=1)
        ax.text(xs[-1], 1.06, "direct annotation", fontsize=7, color=INK2, va="bottom", ha="right")
        ax.set_xscale("log"), ax.set_yscale("log")
        ax.set_xlim(1, nmax * 6)
        xt = [x for x in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000) if x <= nmax]
        ax.set_xticks(xt, [str(x) for x in xt])
        ax.set_ylim(0.1, 100)
        yt = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]
        ax.set_yticks(yt, [f"{t:g}" for t in yt])
        ax.minorticks_off()
        ax.set_xlabel(f"{title} annotated with one self-built gnomAD AF≥1% cache")
        tidy(ax)
        ax.set_title(f"{'AB'[assay == 'wes']}  {title}", loc="left", fontsize=9.5, color=INK, fontweight="bold")
    axes[0].set_ylabel("Time per sample incl. cache build / direct (log)")
    lead_panel(fig.add_subplot(grid[2]), rs)
    fig.text(0.01, 0.01, "A, B: dot = break-even sample count, build / (direct - cached); with a published "
             "cache the build costs nothing.  C: median of 3 samples, published cache (no build).",
             fontsize=7, color=INK2)
    for ext in ("png", "pdf"):
        fig.savefig(RES / f"build_cost_figure.{ext}", dpi=300, facecolor=SURFACE)


def minutes(sec: float) -> str:
    """Compact duration label."""
    m = sec / 60
    if sec < 60:
        return f"{sec:.0f} s"
    if m >= 60:
        return f"{m / 60:.1f} h"
    return f"{m:.0f} min" if m >= 10 else f"{m:.1f} min"


def lead_panel(ax, rs: list[dict]) -> None:
    """Turnaround of one sample: direct against the two cached arms, per annotator and assay."""
    from matplotlib.lines import Line2D

    names = [a[0] for a in ANNOTATORS]
    rows_ = []
    for assay, title in (("wgs", "genome"), ("wes", "exome")):
        for name in names:
            got = {r["cache"]: r for r in rs if r["annotator"] == name and r["assay"] == assay}
            if got:
                rows_.append((f"{SHORT[name]}, {title}", got))
    caches = [("gnomAD AF≥1%", BLUE), ("Cohort", ORANGE)]
    for i, (_, got) in enumerate(rows_):
        y = len(rows_) - 1 - i
        d = next(iter(got.values()))["D"] / 60
        lo = min(got[c]["C"] for c, _ in caches if c in got) / 60
        ax.plot([lo, d], [y, y], color=GRID, lw=2.5, zorder=1, solid_capstyle="round")
        ax.scatter(d, y, s=48, color=SURFACE, edgecolor=INK2, linewidth=1.5, zorder=3)
        ax.annotate(minutes(d * 60), (d, y), xytext=(7, 0), textcoords="offset points", fontsize=7,
                    color=INK2, va="center")
        for c, col in caches:
            if c in got:
                ax.scatter(got[c]["C"] / 60, y, s=48, color=col, edgecolor=SURFACE, linewidth=1.5, zorder=4)
        ax.annotate(minutes(lo * 60), (lo, y), xytext=(-7, 0), textcoords="offset points", fontsize=7,
                    color=INK2, va="center", ha="right")
    ax.set_yticks(range(len(rows_)), [r[0] for r in reversed(rows_)], fontsize=7.8, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.set_xscale("log")
    ax.set_xlim(0.12, 1000)
    xt = [0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500]
    ax.set_xticks(xt, [f"{t:g}" for t in xt])
    ax.minorticks_off()
    ax.set_xlabel("Minutes from called variants to annotated VCF (log)")
    tidy(ax)
    ax.grid(False, axis="y")
    ax.grid(True, axis="x", color=GRID, lw=0.6, zorder=0)
    ax.set_title("C  Turnaround per sample", loc="left", fontsize=9.5, color=INK, fontweight="bold")
    h = [Line2D([], [], ls="", marker="o", ms=7, mfc=SURFACE, mec=INK2, label="direct")] + [
        Line2D([], [], ls="", marker="o", ms=7, color=col, label=lab) for col, lab in
        [(BLUE, "gnomAD AF≥1% cache"), (ORANGE, "cohort cache")]]
    ax.legend(handles=h, loc="lower right", frameon=False, fontsize=7.2, labelcolor=INK)


if __name__ == "__main__":
    rs = rows()
    brows = write_tables(rs, build_costs())
    print((RES / "assay_table.md").read_text())
    print((RES / "build_cost_table.md").read_text())
    if "--no-figure" not in sys.argv:
        assay_figure(rs)
        build_figure(brows, rs)
