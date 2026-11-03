"""The pilot's persistent cache must preserve the existing native semantics."""

import json
import subprocess

from benchmarks.run_pilot import strict_semantic_compare
from tests.test_fused_engine import dataset, native as native_fixture, run_case

native = native_fixture


def test_loaded_cache_two_jobs(native, tmp_path):
    files = dataset(tmp_path / "data", records=40)
    reference, _ = run_case(tmp_path, files, "legacy")
    _, old = run_case(tmp_path, files, "native")
    spec = json.loads((old.work_dir / "native/job.json").read_text())
    jobs = []
    for index in range(2):
        folder = tmp_path / f"loaded-{index}"
        folder.mkdir()
        job = dict(
            spec,
            input=str(files[0]),
            misses=str(folder / "misses.bcf"),
            output=str(folder / "output.bcf"),
            stats=str(folder / "stats.json"),
        )
        # The ordinary native run produced the same annotated miss set; this
        # test exercises cache reuse and output reconstruction, not an annotator.
        jobs.append(job)
    manifest = tmp_path / "jobs.json"
    manifest.write_text(json.dumps(jobs))
    process = subprocess.Popen(
        [native, "resident-batch", "--jobs", str(manifest)],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert json.loads(process.stdout.readline())["event"] == "cache-ready"
    for index in range(2):
        assert json.loads(process.stdout.readline()) == {
            "event": "sample-ready",
            "index": index,
        }
        event = json.loads(process.stdout.readline())
        assert event["event"] == "misses-ready"
        process.stdin.write('{"action":"resume"}\n')
        process.stdin.flush()
        assert json.loads(process.stdout.readline())["event"] == "complete"
        assert strict_semantic_compare(jobs[index]["output"], reference)[
            "semantic_pass"
        ]
    assert json.loads(process.stdout.readline())["event"] == "batch-complete"
    _, err = process.communicate(timeout=15)
    assert process.returncode == 0, err
