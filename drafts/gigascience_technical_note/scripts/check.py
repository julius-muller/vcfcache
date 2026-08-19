#!/usr/bin/env python3
"""Validate manuscript structure, citations, figures, claims, and placeholders."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path

DRAFT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = DRAFT_DIR.parents[1]
PLACEHOLDER = re.compile(r"\[[A-Z][A-Z0-9 ._/,:'()–—-]{3,}\]")
CITATION = re.compile(r"@([A-Za-z0-9_:-]+)")
BIB_KEY = re.compile(r"@[A-Za-z]+\s*\{\s*([^,\s]+)")
WORD = re.compile(r"\b[\w]+(?:[’'-][\w]+)*\b", re.UNICODE)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def words(path: Path) -> int:
    """Count manuscript words after excluding figures and citation markup."""
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"!\[[^]]*\]\([^)]*\)(?:\{[^}]*\})?", "", text)
    text = re.sub(r"\[@[^]]+\]", "", text)
    return len(WORD.findall(text))


def main() -> None:
    """Validate the complete working draft and report publication blockers."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", action="store_true")
    args = parser.parse_args()
    errors: list[str] = []
    warnings: list[str] = []

    required = [
        DRAFT_DIR / "metadata.yaml",
        DRAFT_DIR / "claims.tsv",
        DRAFT_DIR / "references.bib",
        DRAFT_DIR / "figures" / "PROVENANCE.tsv",
    ]
    required.extend(sorted((DRAFT_DIR / "sections").glob("*.md")))
    if len(required) < 15:
        errors.append("Expected complete manuscript section set")
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"Missing or empty required file: {path}")

    abstract_words = words(DRAFT_DIR / "sections" / "01_abstract.md")
    findings_words = words(DRAFT_DIR / "sections" / "03_findings.md")
    methods_words = words(DRAFT_DIR / "sections" / "04_methods.md")
    if not 150 <= abstract_words <= 250:
        errors.append(f"Abstract has {abstract_words} words; expected 150-250")
    if findings_words > 3000:
        errors.append(f"Findings has {findings_words} words; compact ceiling is 3000")
    if methods_words > 2200:
        errors.append(f"Methods has {methods_words} words; compact ceiling is 2200")
    if findings_words < 1100:
        warnings.append(f"Findings is unusually short at {findings_words} words")

    bibliography = (DRAFT_DIR / "references.bib").read_text(encoding="utf-8")
    keys = set(BIB_KEY.findall(bibliography))
    cited: set[str] = set()
    text_files = list((DRAFT_DIR / "sections").glob("*.md"))
    text_files.extend((DRAFT_DIR / "supplement" / "sections").glob("*.md"))
    all_text = "\n".join(path.read_text(encoding="utf-8") for path in text_files)
    cited.update(CITATION.findall(all_text))
    missing_citations = sorted(cited - keys)
    if missing_citations:
        errors.append("Missing bibliography keys: " + ", ".join(missing_citations))

    with (DRAFT_DIR / "claims.tsv").open(encoding="utf-8", newline="") as handle:
        claims = list(csv.DictReader(handle, delimiter="\t"))
    claim_ids = [row["claim_id"] for row in claims]
    if len(claim_ids) != len(set(claim_ids)):
        errors.append("Duplicate claim IDs")
    for row in claims:
        evidence = REPO_ROOT / row["evidence_file"]
        if not evidence.is_file():
            errors.append(f"Claim {row['claim_id']} evidence is missing: {evidence}")
        if row["status"] != "frozen":
            errors.append(f"Claim {row['claim_id']} is not frozen")

    provenance = DRAFT_DIR / "figures" / "PROVENANCE.tsv"
    if provenance.is_file():
        with provenance.open(encoding="utf-8", newline="") as handle:
            figure_rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(figure_rows) < 24:
            errors.append("Expected three-format main and supplementary figure assets")
        for row in figure_rows:
            target = DRAFT_DIR / row["target"]
            if not target.is_file():
                errors.append(f"Figure target is missing: {target}")
            elif sha256(target) != row["sha256"]:
                errors.append(f"Figure checksum changed: {target}")

    placeholders = sorted(set(PLACEHOLDER.findall(all_text + "\n" + (DRAFT_DIR / "metadata.yaml").read_text(encoding="utf-8"))))
    if placeholders:
        message = f"{len(placeholders)} author-controlled placeholders remain"
        if args.release:
            errors.append(message + ": " + "; ".join(placeholders))
        else:
            warnings.append(message + " (expected in a working draft)")

    print(f"Abstract words: {abstract_words}")
    print(f"Findings words: {findings_words}")
    print(f"Methods words: {methods_words}")
    print(f"Claims: {len(claims)}")
    print(f"Citations: {len(cited)} used / {len(keys)} defined")
    for warning in warnings:
        print(f"WARNING: {warning}")
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        raise SystemExit(1)
    print("Manuscript checks passed")


if __name__ == "__main__":
    main()
