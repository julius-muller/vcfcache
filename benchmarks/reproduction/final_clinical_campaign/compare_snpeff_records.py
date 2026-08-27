#!/usr/bin/env python3
"""Compare SnpEff outputs and distinguish ordering from content differences."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def rows(path: Path) -> dict[tuple[str, str, str, str], tuple[str, str, str]]:
    fmt = r"%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t%INFO/CLNSIG\t%INFO/AF\n"
    proc = subprocess.Popen(
        ["bcftools", "query", "--allow-undef-tags", "-f", fmt, str(path)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    result: dict[tuple[str, str, str, str], tuple[str, str, str]] = {}
    for line in proc.stdout:
        chrom, pos, ref, alt, ann, clin_sig, af = line.rstrip("\n").split("\t")
        key = chrom, pos, ref, alt
        if key in result:
            raise ValueError(f"duplicate variant key: {key}")
        result[key] = ann, clin_sig, af
    if proc.wait() != 0:
        raise RuntimeError(f"bcftools query failed for {path}")
    return result


def canonical(value: str) -> tuple[str, ...]:
    return tuple(sorted(value.split(",")))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("direct", type=Path)
    parser.add_argument("cached", type=Path)
    parser.add_argument("--limit", type=int, default=100)
    args = parser.parse_args()

    direct = rows(args.direct)
    cached = rows(args.cached)
    direct_only = sorted(direct.keys() - cached.keys())
    cached_only = sorted(cached.keys() - direct.keys())
    raw_differences = 0
    content_differences = 0
    examples: list[str] = []
    for key in sorted(direct.keys() & cached.keys()):
        d_ann, d_clin, d_af = direct[key]
        c_ann, c_clin, c_af = cached[key]
        if (d_ann, d_clin, d_af) == (c_ann, c_clin, c_af):
            continue
        raw_differences += 1
        ordering_only = (
            canonical(d_ann) == canonical(c_ann)
            and canonical(d_clin) == canonical(c_clin)
            and canonical(d_af) == canonical(c_af)
        )
        if not ordering_only:
            content_differences += 1
        if len(examples) < args.limit:
            examples.append(
                "\t".join(
                    [
                        *key,
                        "ORDER_ONLY" if ordering_only else "CONTENT",
                        d_ann,
                        c_ann,
                        d_clin,
                        c_clin,
                        d_af,
                        c_af,
                    ]
                )
            )

    print(f"direct_records\t{len(direct)}")
    print(f"cached_records\t{len(cached)}")
    print(f"direct_only\t{len(direct_only)}")
    print(f"cached_only\t{len(cached_only)}")
    print(f"raw_value_differences\t{raw_differences}")
    print(
        "content_differences_after_list_order_canonicalization"
        f"\t{content_differences}"
    )
    print(
        "#CHROM\tPOS\tREF\tALT\tCLASS\tDIRECT_ANN\tCACHED_ANN\t"
        "DIRECT_CLNSIG\tCACHED_CLNSIG\tDIRECT_AF\tCACHED_AF"
    )
    print("\n".join(examples))


if __name__ == "__main__":
    main()
