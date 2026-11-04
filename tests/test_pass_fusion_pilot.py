"""Guard branch-specific direct arguments and durable review behaviour."""

import json

import pytest

from benchmarks.pass_fusion_pilot.run import MODES, collect, command_for
from benchmarks.run_pilot import sha256sum, write_json_atomic


@pytest.mark.parametrize("mode", MODES)
def test_cli_control_engine(mode, tmp_path):
    """This Python branch rejects adaptive routing combined with --uncached."""
    command = command_for(mode, *([tmp_path / "file"] * 5))
    assert command[command.index("--engine") + 1] == (
        "legacy" if mode == "direct" else "adaptive"
    )
    assert ("--uncached" in command) == (mode == "direct")
    assert ("--engine-workers" in command) == (mode not in {"direct", "control-auto"})


def test_collect_preserves_case_and_nonterminal_states(tmp_path):
    """A runner state already contains a case; review must not duplicate kwargs."""
    cases = [dict(id="finished", node="node-a"), dict(id="pending", node="node-b")]
    write_json_atomic(tmp_path / "manifest.json", dict(cases=cases))
    (tmp_path / "READY").write_text(sha256sum(tmp_path / "manifest.json") + "\n")
    path = tmp_path / "cases/finished"
    path.mkdir(parents=True)
    write_json_atomic(
        path / "state.json", dict(case=cases[0], state="complete", rows=[])
    )
    collect(tmp_path)
    result = json.loads((tmp_path / "REVIEW.json").read_text())
    assert result["complete"] == 1 and result["not_terminal"] == 1
    assert result["cases"][0]["case"] == cases[0]
