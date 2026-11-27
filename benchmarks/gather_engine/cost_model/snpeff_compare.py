"""Compare two SnpEff outputs record by record, in file order.

Each record is identical, identical up to the order of INFO tags, differs only in the warnings
field of ANN (SnpEff 5.4c adds WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS depending on which records
it processed before; benchmarks/docs/SNPEFF_PROCESSING_ORDER_BUG.md), or differs otherwise.
Prints the counts, the warnings involved and examples of any other difference.
Usage: snpeff_compare.py direct.bcf cached.bcf
"""

import collections
import itertools
import subprocess
import sys

WARN = 15  # "ERRORS / WARNINGS / INFO", the last of the 16 ANN fields


def records(path):
    """Yield the tab-split body lines of a BCF."""
    with subprocess.Popen(["bcftools", "view", "-H", path], stdout=subprocess.PIPE, text=True) as p:
        for line in p.stdout:
            yield line.rstrip("\n").split("\t")
        if p.wait():
            raise SystemExit(f"{path}: bcftools failed")


def info(rec):
    """Return the INFO column as a dict."""
    return dict(x.split("=", 1) if "=" in x else (x, "") for x in rec[7].split(";"))


def ann_diff(a, b):
    """Return the set of warning tokens that differ if the ANN values differ only there, else None."""
    ea, eb = a.split(","), b.split(",")
    if len(ea) != len(eb):
        return None
    tokens = set()
    for x, y in zip(ea, eb, strict=True):
        u, v = x.split("|"), y.split("|")
        if len(u) != len(v) or u[:WARN] != v[:WARN] or u[WARN + 1 :] != v[WARN + 1 :]:
            return None
        tokens |= set(u[WARN].split("&")) ^ set(v[WARN].split("&"))
    return tokens


def main():
    """Compare the two files and print a report."""
    direct, cached = sys.argv[1], sys.argv[2]
    counts, warnings, examples = collections.Counter(), collections.Counter(), []
    end = object()
    for a, b in itertools.zip_longest(records(direct), records(cached), fillvalue=end):
        if a is end or b is end:
            raise SystemExit("record counts differ")
        if a[:5] != b[:5]:
            raise SystemExit(f"records out of step: {a[:5]} vs {b[:5]}")
        counts["records"] += 1
        if a == b:
            counts["identical"] += 1
            continue
        ia, ib = info(a), info(b)
        if a[:7] == b[:7] and ia == ib:
            counts["identical (INFO tag order only)"] += 1
            continue
        other = {k for k in ia.keys() | ib.keys() if k != "ANN" and ia.get(k) != ib.get(k)}
        if a[:7] != b[:7]:
            other.add("<columns 1-7>")
        tokens = ann_diff(ia.get("ANN", ""), ib.get("ANN", "")) if ia.get("ANN") != ib.get("ANN") else set()
        if not other and tokens:
            counts["ANN warnings only"] += 1
            warnings.update(tokens)
        else:
            counts["other"] += 1
            if len(examples) < 10:
                fields = other | ({"ANN"} if tokens is None else set())
                examples.append(f"{':'.join(a[:2])} {a[3]}>{a[4]}: {sorted(fields)}")
    for k, v in counts.items():
        print(f"{k}\t{v}")
    if warnings:
        print("ANN warnings\t" + ",".join(f"{w}:{n}" for w, n in warnings.most_common(10)))
    for e in examples:
        print("example\t" + e)


if __name__ == "__main__":
    main()
