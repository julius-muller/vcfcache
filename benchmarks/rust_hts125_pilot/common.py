import json
import shlex
import subprocess
import tempfile
from pathlib import Path

from benchmarks.run_pilot import sha256sum


def mapped(path):
    """Controller and worker spellings of the shared mount."""
    return str(path).replace("/mnt/data/slurm-results/", "/results/", 1)


def controller(path):
    """Resolve a worker-visible path on the controller."""
    return Path(str(path).replace("/results/", "/mnt/data/slurm-results/", 1))


def definition(root):
    """Require the frozen definition, never an edited live job recipe."""
    assert sha256sum(root / "manifest.json") == (root / "READY").read_text().strip()
    return json.loads((root / "manifest.json").read_text())


def verify(path, digest):
    """Verify immutable inputs/cache payloads without editing them."""
    path = Path(path)
    assert path.exists() and sha256sum(path) == digest, f"Asset changed: {path}"
    return path


def wrapper(directory, actual, latest):
    """Thread every supported stage, including commands inside frozen recipes."""
    directory.mkdir()
    destination = directory / "bcftools"
    script = f"""#!/bin/bash
set -euo pipefail
actual={shlex.quote(str(actual))}
subcommand=${{1-}}
case "$subcommand" in
  annotate|view|filter|concat|index{'|sort' if latest else ''})
    shift
    for arg in "$@"; do
      case "$arg" in --threads|--threads=*|-@*) exec "$actual" "$subcommand" "$@";; esac
    done
    exec "$actual" "$subcommand" "$@" --threads 8;;
esac
exec "$actual" "$@"
"""
    destination.write_text(script)
    destination.chmod(0o755)
    return destination


def prepared_input(task, scratch):
    """Assay extraction/BCF encoding only; never split or remove any variants."""
    source = verify(task["input"]["path"], task["input"]["sha256"])
    output = scratch / (task["assay"] + ".bcf")
    command = [
        "/usr/local/bin/bcftools",
        "view",
        str(source),
        "-Ob1",
        "-W",
        "--threads",
        "8",
        "-o",
        str(output),
    ]
    if task.get("preprocess", {}).get("bed"):
        bed = task["preprocess"]["bed"]
        command += ["-R", str(verify(bed["path"], bed["sha256"]))]
    subprocess.run(command, check=True)
    return output


def preflight_prepared(bcftools, source, tag):
    """Verify the experimental trusted-input assertion outside all timings."""
    header = subprocess.check_output([bcftools, "view", "-h", str(source)], text=True)
    has_tag = f"##INFO=<ID={tag}," in header
    fmt = "%CHROM\t%POS\t%REF\t%ALT" + (f"\t%INFO/{tag}" if has_tag else "") + "\n"
    locus, seen, records = None, set(), 0
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(
            [bcftools, "query", "-f", fmt, str(source)],
            stdout=subprocess.PIPE,
            stderr=errors,
            text=True,
        )
        assert process.stdout is not None
        try:
            with process.stdout:
                for line in process.stdout:
                    fields = line.rstrip("\n").split("\t")
                    contig, pos, ref, alt = fields[:4]
                    if "," in alt or alt in {"*", ".", ""}:
                        raise ValueError(
                            "Pilot requires already-biallelic, star-free input; do not silently transform it"
                        )
                    if has_tag and fields[4] not in {"", "."}:
                        raise ValueError(
                            "Existing input annotation could masquerade as a cache hit"
                        )
                    if locus != (contig, pos):
                        locus, seen = (contig, pos), set()
                    if (ref, alt) in seen:
                        raise ValueError(
                            "Duplicate variant keys are outside this prepared-input pilot"
                        )
                    seen.add((ref, alt))
                    records += 1
            if process.wait():
                errors.seek(0)
                raise RuntimeError(errors.read().decode(errors="replace"))
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait()
    return dict(
        records=records,
        biallelic=True,
        star_free=True,
        unique_keys=True,
        unannotated=True,
    )
