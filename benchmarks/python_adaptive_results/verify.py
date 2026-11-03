"""Offline checks of the archived evidence, exported measurements and final figures."""

import argparse
import csv
import hashlib
import json
import tarfile
from collections import Counter
from pathlib import Path

from benchmarks.python_adaptive_results.collect import flatten


def verify(root):
    """Check all 600 exported rows against their original evidence, offline."""
    snapshot = json.loads((root / "source/SNAPSHOT.json").read_text())
    archive = root / "source/terminal-evidence.tar.gz"
    assert (
        hashlib.sha256(archive.read_bytes()).hexdigest() == snapshot["archive_sha256"]
    )
    # Inspect the archive before any extraction: only simple regular files and
    # directories are accepted. Extraction never touches the working snapshot.
    import tempfile

    with tempfile.TemporaryDirectory(prefix="vcfcache-verify-") as directory:
        with tarfile.open(archive) as stream:
            for member in stream.getmembers():
                assert (
                    not member.name.startswith("/")
                    and ".." not in Path(member.name).parts
                )
                assert member.isfile() or member.isdir()
            stream.extractall(directory, filter="data")
        campaign = Path(directory) / "campaign"
        manifest = campaign / "manifest.json"
        assert (
            hashlib.sha256(manifest.read_bytes()).hexdigest()
            == snapshot["manifest_sha256"]
        )
        data = json.loads(manifest.read_text())
        latest = {}
        for p in sorted((campaign / "attempts").glob("*/*/state.json")):
            latest[p.parent.name] = p.parent
        rows = [flatten(t, latest[t["id"]], manifest) for t in data["tasks"]]
        with (root / "tables/cells.tsv").open() as stream:
            exported = list(csv.DictReader(stream, delimiter="\t"))
        # Compare every exported field, including identity, scope and timings.
        assert len(rows) == len(exported) == 600
        for row, saved in zip(rows, exported, strict=True):
            assert {k: "" if v is None else str(v) for k, v in row.items()} == saved
        assert Counter(r["state"] for r in rows) == {
            "complete": 572,
            "complete_with_direct_difference": 24,
            "failed": 4,
        }
        with (root / "tables/plot_views.tsv").open() as stream:
            views = list(csv.DictReader(stream, delimiter="\t"))
        canonical = {r["condition_id"]: r for r in exported}
        assert len(views) == 612
        for v in views:
            c = canonical[v["measurement_id"]]
            for field in [
                "direct_seconds",
                "cached_seconds",
                "input_sha256",
                "cache_sha256",
                "job_id",
            ]:
                assert v[field] == c[field]
        figures = root / "figures-v2"
        for stem in [
            "manuscript/Figure1_workflow",
            "manuscript/Figure2_empirical_performance",
            "manuscript/Figure3_WGS_cache_strategies",
            "supplements/FigureS1_pipeline_complexity",
            "supplements/FigureS2_assays_by_assembly",
            "supplements/FigureS3_resources",
            "repository/WGS_time_saved",
        ]:
            for suffix in [".pdf", ".svg", ".png"]:
                p = figures / (stem + suffix)
                assert p.is_file() and p.stat().st_size > 1000
        for panel in [
            "Figure2A_assay_speedup",
            "Figure2B_WGS_time_saved",
            "Figure2C_WGS_throughput",
            "Figure3A_WGS_cache_hits",
            "Figure3B_WGS_speedup_distributions",
        ]:
            assert (figures / "manuscript/panels" / (panel + ".pdf")).is_file()
    print(
        "PASS: archive SHA256; all 600 rows re-audited; 612 views; no alias timing duplication; final figures and panels present"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    verify(parser.parse_args().snapshot)
