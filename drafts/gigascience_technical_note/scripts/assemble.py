#!/usr/bin/env python3
"""Assemble independently editable manuscript sections into draft files."""

from __future__ import annotations

from pathlib import Path

DRAFT_DIR = Path(__file__).resolve().parents[1]
BUILD_DIR = DRAFT_DIR / "build"


def assemble(
    source_dir: Path,
    output: Path,
    metadata: Path,
    *,
    exclude: frozenset[str] = frozenset(),
) -> None:
    """Combine ordered section files with metadata into one Markdown draft."""
    parts: list[str] = []
    parts.append(metadata.read_text(encoding="utf-8").strip())
    for path in sorted(source_dir.glob("*.md")):
        if path.name in exclude:
            continue
        parts.append(path.read_text(encoding="utf-8").strip())
    output.write_text("\n\n".join(parts) + "\n", encoding="utf-8")


def main() -> None:
    """Assemble the main manuscript and supplementary material."""
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    assemble(
        DRAFT_DIR / "sections",
        BUILD_DIR / "vcfcache_gigascience_technical_note.md",
        DRAFT_DIR / "metadata.yaml",
        exclude=frozenset({"11_figure_legends.md"}),
    )
    assemble(
        DRAFT_DIR / "supplement" / "sections",
        BUILD_DIR / "vcfcache_supplementary_material.md",
        DRAFT_DIR / "supplement" / "metadata.yaml",
    )


if __name__ == "__main__":
    main()
