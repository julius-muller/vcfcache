"""Compare two Ensembl VEP outputs record by record, in file order.

Each record is identical, differs only in HGNC_ID (the known VEP bug, Ensembl/ensembl-vep#1959),
differs only in the order of a list, or differs otherwise. VEP lists are '&'-separated, and CLIN_SIG
joins several ClinVar submissions with ';' (written %3B) in Perl hash order, which changes from run
to run even for one variant run alone. Prints the counts, the genes behind the HGNC_ID differences
and examples of any other difference.
Usage: vep116_compare.py direct.bcf cached.bcf
"""

import collections
import itertools
import subprocess
import sys


def csq_format(path):
    """Return the field names of the CSQ tag."""
    head = subprocess.run(["bcftools", "view", "-h", path], capture_output=True, text=True, check=True)
    for line in head.stdout.splitlines():
        if line.startswith("##INFO=<ID=CSQ"):
            return line.split("Format: ")[1].rstrip('">').split("|")
    raise SystemExit(f"{path}: no CSQ header")


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


def unordered(value):
    """Return a CSQ value with the order of its %3B groups and &-lists removed."""
    return sorted(tuple(sorted(g.split("&"))) for g in value.split("%3B"))


def csq_diff(a, b, fmt):
    """Return the CSQ fields that differ, and those of them that differ in more than list order."""
    ea, eb = a.split(","), b.split(",")
    if len(ea) != len(eb):
        return {"<CSQ entry count>"}, {"<CSQ entry count>"}
    fields, real = set(), set()
    for x, y in zip(ea, eb, strict=True):
        for k, u, v in zip(fmt, x.split("|"), y.split("|"), strict=True):
            if u != v:
                fields.add(k)
                if unordered(u) != unordered(v):
                    real.add(k)
    return fields, real


def main():
    """Compare the two files and print a report."""
    direct, cached = sys.argv[1], sys.argv[2]
    fmt = csq_format(direct)
    if csq_format(cached) != fmt:
        raise SystemExit("CSQ formats differ")
    sym = fmt.index("SYMBOL")
    counts, genes, examples = collections.Counter(), collections.Counter(), []
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
        other = {k for k in ia.keys() | ib.keys() if k != "CSQ" and ia.get(k) != ib.get(k)}
        fields, real = csq_diff(ia.get("CSQ", ""), ib.get("CSQ", ""), fmt)
        if a[:7] != b[:7]:
            other.add("<columns 1-7>")
        if not other and real <= {"HGNC_ID"} and "HGNC_ID" in fields:
            counts["HGNC_ID only" if fields == {"HGNC_ID"} else "HGNC_ID and list order only"] += 1
            for x in ia["CSQ"].split(","):
                genes[x.split("|")[sym]] += 1
        elif not other and not real:
            counts["list order only: " + ",".join(sorted(fields))] += 1
        else:
            counts["other"] += 1
            if len(examples) < 10:
                examples.append(f"{':'.join(a[:2])} {a[3]}>{a[4]}: {sorted(other | fields)}")
    for k, v in counts.items():
        print(f"{k}\t{v}")
    if genes:
        print("HGNC_ID genes\t" + ",".join(g or "-" for g, _ in genes.most_common(30)))
    for e in examples:
        print("example\t" + e)


if __name__ == "__main__":
    main()
