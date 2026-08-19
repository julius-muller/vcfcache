#!/usr/bin/env python3
"""Build editable DOCX and review PDF files from the Markdown sources."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from assemble import BUILD_DIR, DRAFT_DIR
from assemble import main as assemble

DOCUMENTS = (
    ("vcfcache_gigascience_technical_note", "metadata.yaml"),
    ("vcfcache_supplementary_material", "supplement/metadata.yaml"),
)


def require_program(name: str) -> str:
    """Return an executable path or stop with a concise build error."""
    path = shutil.which(name)
    if path is None:
        raise SystemExit(f"Required program not found: {name}")
    return path


def run(command: list[str], env: dict[str, str] | None = None) -> None:
    """Run one checked manuscript-build command from the draft root."""
    subprocess.run(command, cwd=DRAFT_DIR, env=env, check=True)


def file_sha256(path: Path) -> str:
    """Calculate the SHA-256 digest of one build artifact."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


WORD_NAMESPACE = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
WORD = f"{{{WORD_NAMESPACE}}}"


def polish_docx(path: Path) -> None:
    """Apply journal-like body typography while retaining an editable DOCX."""
    with zipfile.ZipFile(path) as archive:
        members = {name: archive.read(name) for name in archive.namelist()}

    root = ET.fromstring(members["word/styles.xml"])
    for style in root.findall(f"{WORD}style"):
        style_id = style.get(f"{WORD}styleId")
        if style_id in {"Normal", "BodyText", "FirstParagraph", "Bibliography"}:
            paragraph_properties = style.find(f"{WORD}pPr")
            if paragraph_properties is None:
                paragraph_properties = ET.Element(f"{WORD}pPr")
                style.insert(0, paragraph_properties)
            justification = paragraph_properties.find(f"{WORD}jc")
            if justification is None:
                justification = ET.SubElement(paragraph_properties, f"{WORD}jc")
            justification.set(f"{WORD}val", "both")

        if style_id in {"Heading1", "Heading2", "Heading3"}:
            run_properties = style.find(f"{WORD}rPr")
            if run_properties is None:
                run_properties = ET.SubElement(style, f"{WORD}rPr")
            colour = run_properties.find(f"{WORD}color")
            if colour is None:
                colour = ET.SubElement(run_properties, f"{WORD}color")
            colour.attrib.clear()
            colour.set(f"{WORD}val", "000000")

    ET.register_namespace("w", WORD_NAMESPACE)
    members["word/styles.xml"] = ET.tostring(
        root, encoding="utf-8", xml_declaration=True
    )

    temporary = path.with_suffix(".polished.docx")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in members.items():
            archive.writestr(name, content)
    temporary.replace(path)


def main() -> None:
    """Assemble and render the manuscript and supplementary document."""
    assemble()
    pandoc = require_program("pandoc")
    libreoffice = require_program("libreoffice")

    for stem, metadata in DOCUMENTS:
        source = BUILD_DIR / f"{stem}.md"
        docx = BUILD_DIR / f"{stem}.docx"
        run(
            [
                pandoc,
                "--standalone",
                "--citeproc",
                "--csl=styles/gigascience.csl",
                f"--metadata-file={metadata}",
                "--bibliography=references.bib",
                "--resource-path=.",
                str(source.relative_to(DRAFT_DIR)),
                "--output",
                str(docx.relative_to(DRAFT_DIR)),
            ]
        )
        polish_docx(docx)

        latex = BUILD_DIR / f"{stem}.tex"
        run(
            [
                pandoc,
                "--standalone",
                "--citeproc",
                "--csl=styles/gigascience.csl",
                f"--metadata-file={metadata}",
                "--bibliography=references.bib",
                "--resource-path=.",
                str(source.relative_to(DRAFT_DIR)),
                "--output",
                str(latex.relative_to(DRAFT_DIR)),
            ]
        )

    with tempfile.TemporaryDirectory(prefix="vcfcache-manuscript-lo-") as profile:
        profile_uri = Path(profile).resolve().as_uri()
        environment = os.environ.copy()
        docx_files = [str(BUILD_DIR / f"{stem}.docx") for stem, _metadata in DOCUMENTS]
        run(
            [
                libreoffice,
                "--headless",
                f"-env:UserInstallation={profile_uri}",
                "--convert-to",
                "pdf",
                "--outdir",
                str(BUILD_DIR),
                *docx_files,
            ],
            env=environment,
        )

    artifacts = sorted(
        path
        for path in BUILD_DIR.iterdir()
        if path.suffix in {".md", ".tex", ".docx", ".pdf"}
    )
    checksums = "\n".join(
        f"{file_sha256(path)}  {path.name}" for path in artifacts
    )
    (BUILD_DIR / "CHECKSUMS.sha256").write_text(checksums + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
