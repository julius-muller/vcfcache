"""Thread propagation and durable partial-review regression checks."""

import json
import subprocess

from benchmarks.input_preserving_pilot.run import INTEL, collect, wrapper
from benchmarks.run_pilot import sha256sum


def test_thread_wrapper_covers_recipe_stages(tmp_path):
    """Supported stages gain threads; explicit settings and query stay intact."""
    actual = tmp_path / "actual"
    actual.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\n')
    actual.chmod(0o755)
    latest = wrapper(tmp_path / "latest", actual, True)
    old = wrapper(tmp_path / "old", actual, False)
    for command in ("view", "annotate", "filter", "concat", "index", "sort"):
        args = subprocess.check_output(
            [str(latest), command, "input"], text=True
        ).splitlines()
        assert args[-2:] == ["--threads", "8"]
    args = subprocess.check_output(
        [str(latest), "view", "--threads", "2", "input"], text=True
    ).splitlines()
    assert args == ["view", "--threads", "2", "input"]
    assert subprocess.check_output(
        [str(latest), "query", "input"], text=True
    ).splitlines() == ["query", "input"]
    assert subprocess.check_output(
        [str(old), "sort", "input"], text=True
    ).splitlines() == ["sort", "input"]
    assert "sl-w5" not in INTEL and len(INTEL) == 5


def test_collector_preserves_nonterminal_evidence(tmp_path):
    """A review is not an execution or permission to submit more work."""
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(dict(cases=[dict(id="test", node="sl-w1")], deployment="unused"))
    )
    (tmp_path / "READY").write_text(sha256sum(manifest))
    folder = tmp_path / "cases/test"
    folder.mkdir(parents=True)
    original = json.dumps(
        dict(state="running", rows=[dict(mode="new-direct", seconds=1.0, records=1)])
    )
    (folder / "state.json").write_text(original)
    collect(tmp_path)
    assert (folder / "state.json").read_text() == original
    assert "running" in (tmp_path / "RESULTS.md").read_text()
    assert not (tmp_path / "terminal-evidence.tar.gz").exists()
