"""Resolve historical worker identities; read-only, never submit benchmarks.

Only --refresh uses SSH. Rebuilds use the frozen compact provenance snapshot.
Unknown placement is not silently classified as Intel.
"""

import csv
import hashlib
import io
import json
import re
import subprocess
import tarfile
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
V1 = ROOT.parent / "decision_20261109"
PY = "367100833"
RUST = "e4d278ce1"
PROOFS = {}


def digest(raw):
    """Hash complete source bytes for audit provenance."""
    return hashlib.sha256(raw).hexdigest()


def host(value):
    """Worker identity in a resource cgroup, not a logical group such as g00."""
    found = set(re.findall(r"(?<![a-zA-Z0-9])sl-w[1-6](?![0-9])", str(value)))
    return next(iter(found)) if len(found) == 1 else ""


@lru_cache(None)
def git(path, commit=RUST):
    """Read and hash an immutable Git object without touching worktrees."""
    raw = subprocess.check_output(["git", "show", commit + ":" + path], cwd=REPO)
    PROOFS[commit + ":" + path] = digest(raw)
    return raw


def table(path, commit=RUST):
    """Parse a pinned TSV export."""
    return list(csv.DictReader(io.StringIO(git(path, commit).decode()), delimiter="\t"))


def json_git(path):
    """Read pinned resource metadata, retaining missing proofs as unknown."""
    try:
        return json.loads(git(path))
    except subprocess.CalledProcessError:
        return {}


@lru_cache(None)
def archive(path, commit=RUST):
    """Read JSON members in memory; never extract potentially unsafe paths."""
    with tarfile.open(fileobj=io.BytesIO(git(path, commit))) as stream:
        return {
            m.name: json.load(stream.extractfile(m))
            for m in stream.getmembers()
            if m.isfile()
            and m.name.endswith(("resources.json", "state.json", "manifest.json"))
        }


def key(row):
    """Identify the same observation before and after numeric serialization."""
    return "|".join(
        [
            row["config_id"],
            row["condition_id"],
            row["sample"],
            str(int(float(row["repetition"]))),
        ]
    )


def resource_host(value):
    """Extract worker identity from a condition's recorded cgroup path."""
    return host(value.get("cgroup_after", {}).get("path", ""))


def refresh():
    """Freeze worker proofs for all v1 observations, including direct references."""
    configs = {r["config_id"]: r for r in read_tsv(V1 / "configurations.tsv")}
    rows = read_tsv(V1 / "measurements.tsv")
    rust_rows = table(
        "benchmarks/rust_adaptive_results/final-20261103/tables/plot_views.tsv", PY
    )
    rust = {r["measurement_id"]: r for r in rust_rows}
    matrix = table(
        "benchmarks/branch_eval/integration-review-20261023/cells_backend_audit.tsv"
    )
    wes = archive("benchmarks/wes_pilot/results-20261103/terminal-evidence.tar.gz")
    decoder = archive(
        "benchmarks/wes_pilot/decoder-results-20261103/terminal-evidence.tar.gz"
    )
    manifests = archive(
        "benchmarks/python_adaptive_results/final-20261103/source/terminal-evidence.tar.gz",
        PY,
    )
    tasks = {t["id"]: t for t in manifests["campaign/manifest.json"]["tasks"]}
    profile = table(
        "benchmarks/engine_profile/consolidated-20261027/canonical_conditions.tsv"
    )
    result, requests = {}, {}

    def wes_proof(campaign, cell, mode, bundle=wes):
        prefix = "source/" + (campaign + "/" if campaign else "") + "cases/" + cell
        resource = bundle.get(prefix + "/" + mode + "/resources.json", {})
        state = bundle.get(prefix + "/state.json", {})
        direct = bundle.get(prefix + "/direct/resources.json", {})
        if not direct:
            reuse = state.get("reused_direct", {}).get("path", "")
            match = re.search(r"/wes-pilot/(.+)$", reuse)
            if match:
                direct = wes.get("source/" + match[1] + "/direct/resources.json", {})
        return (
            resource_host(resource),
            resource_host(direct),
            prefix,
            direct.get("wall_seconds"),
        )

    for r in rows:
        group = configs[r["config_id"]]["plot_group"]
        p = dict(
            node=host(r["node"]),
            direct_node="",
            cached_proof=r["source_id"],
            direct_proof="",
            direct_seconds=None,
        )
        if group in {"cohort-primary", "cohort-pgp", "wgs-strategies", "delay"}:
            original = rust[r["condition_id"]]
            p["direct_node"] = host(original.get("direct_step", ""))
            p["direct_proof"] = "Rust plot_views.tsv direct_step: " + r["condition_id"]
            p["direct_seconds"] = float(original["direct_seconds"])
            if not p["direct_node"]:
                task = tasks[r["condition_id"]]
                baseline = task["baseline"]
                artifact = baseline.get("direct", {}).get("path", "")
                if artifact and group != "delay":
                    p["remote_reference"] = artifact
                    requests[artifact] = float(r["direct_s"])
        elif group in {"wes-cohort", "wes-workers"}:
            campaign = r["source_location"].split(":")[0]
            cell, mode = r["condition_id"].split(":")
            n, d, proof, seconds = wes_proof(campaign, cell, mode)
            p.update(
                node=n,
                direct_node=d,
                cached_proof=proof,
                direct_proof=proof + " state/direct",
                direct_seconds=seconds,
            )
        elif group == "decoder":
            cell, mode = r["condition_id"].split(":")
            n, d, proof, seconds = wes_proof("", cell, mode, decoder)
            p.update(
                node=n,
                direct_node=d,
                cached_proof="decoder archive:" + proof,
                direct_proof=proof + " reused_direct",
                direct_seconds=seconds,
            )
        elif group == "matrix":
            job, scale, arm = r["condition_id"].split(":")
            original = next(
                x
                for x in matrix
                if x["job"] == job and x["scale"] == scale and x["requested_arm"] == arm
            )
            candidates = [
                x
                for x in matrix
                if x["primary"].lower() == "true"
                and x["requested_arm"] == "direct"
                and x["recipe"] == original["recipe"]
                and x["scale"] == scale
            ]
            direct = next(
                (x for x in candidates if x["job"] == job),
                candidates[0] if len(candidates) == 1 else None,
            )

            def conditions(x):
                path = (
                    "benchmarks/branch_eval/results-20261023/"
                    + x["comparison"].split("/")[0]
                    + "/CONDITIONS.txt"
                )
                return (
                    host(
                        next(
                            (
                                s
                                for s in git(path).decode().splitlines()
                                if s.startswith("node=")
                            ),
                            "",
                        )
                    ),
                    path,
                )

            p["node"], p["cached_proof"] = conditions(original)
            if direct:
                p["direct_node"], p["direct_proof"] = conditions(direct)
                p["direct_seconds"] = float(direct["wall_seconds"])
        elif group == "profile":
            campaign, cell = r["condition_id"].split(":")
            original = next(
                x for x in profile if x["campaign"] == campaign and x["cell"] == cell
            )
            prefix = "benchmarks/engine_profile/consolidated-20261027/source/"
            path = (
                prefix
                + campaign
                + "/attempts/"
                + original["attempt"]
                + "/"
                + cell
                + "/resources.json"
            )
            p.update(node=resource_host(json_git(path)), cached_proof=path)
            candidates = [
                x
                for x in profile
                if x["campaign"] == campaign
                and x["group"] == original["group"]
                and x["assay"] == original["assay"]
                and x["recipe"] == original["recipe"]
                and x["role"] == "direct"
                and x["state"] == "complete"
            ]
            if len(candidates) == 1:
                d = candidates[0]
                path = (
                    prefix
                    + campaign
                    + "/attempts/"
                    + d["attempt"]
                    + "/"
                    + d["cell"]
                    + "/resources.json"
                )
                p.update(
                    direct_node=resource_host(json_git(path)),
                    direct_proof=path,
                    direct_seconds=float(d["wall_seconds"]),
                )
        elif group in {"fusion", "main125", "python125", "rust125"}:
            # These blocks execute direct and cached on the recorded same worker.
            p.update(
                direct_node=p["node"],
                direct_proof=r["source_id"] + " same-case block",
                direct_seconds=float(r["direct_s"]) if r["direct_s"] else None,
            )
        result[key(r)] = p

    # Only tiny JSON metadata is read from the controller. No benchmark commands,
    # scheduling changes, or writes are performed remotely.
    script = """
import json,pathlib,hashlib
requests = REQUESTS
out={}
for artifact in requests:
 root=pathlib.Path(artifact.replace('/results/','/mnt/data/slurm-results/',1)).parent
 f=root/'metrics.json'
 try:
  raw=f.read_bytes(); m=json.loads(raw)
  out[artifact]=dict(path=str(f),sha256=hashlib.sha256(raw).hexdigest(),wall_seconds=m.get('wall_seconds'),cgroup=m.get('cgroup_v2',{}).get('after',{}).get('path',''))
 except (OSError,ValueError) as e: out[artifact]=dict(error=str(e))
print(json.dumps(out))
""".replace("REQUESTS", repr(list(requests)))
    remote = json.loads(
        subprocess.check_output(
            [
                "ssh",
                "-o",
                "BatchMode=yes",
                "-o",
                "ConnectTimeout=10",
                "ubuntu@10.133.241.128",
                "python3",
                "-",
            ],
            input=script.encode(),
        )
    )
    for p in result.values():
        if "remote_reference" in p:
            ref = remote[p["remote_reference"]]
            p.update(
                direct_node=host(ref.get("cgroup", "")),
                direct_proof=ref.get("path", ""),
                direct_seconds=ref.get("wall_seconds"),
            )
    snapshot = dict(
        rows=result,
        remote=remote,
        source_sha256=PROOFS,
        policy="No sl-w5 cached run or denominator; unknown placement excluded; paired cohort intersection",
    )
    (ROOT / "node_provenance.json").write_text(json.dumps(snapshot, indent=2) + "\n")
    return snapshot


def read_tsv(path):
    """Read a local frozen TSV without mutating its contents."""
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def load(refresh_metadata=False):
    """Prefer frozen local provenance; SSH is an explicit metadata refresh."""
    path = ROOT / "node_provenance.json"
    return (
        refresh()
        if refresh_metadata or not path.exists()
        else json.loads(path.read_text())
    )
