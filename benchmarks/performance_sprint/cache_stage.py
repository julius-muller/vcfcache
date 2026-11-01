"""Byte-identical fallback for drifted worker-local cache replicas.

Keep the original local/shared storage class. Canonical payload transport and
one-time node installation are preparation, outside timed annotation cells.
Original caches are never overwritten.
"""

import fcntl
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

from benchmarks.run_pilot import sha256sum, write_json_atomic


def snapshot_metadata(asset, files, output):
    """Freeze canonical cache provenance, not another node's independent build log."""
    path = Path(asset["path"])
    if not str(path).startswith("/mnt/data/") or "provenance" not in path.name:
        return asset
    destination = output / "cache-provenance" / asset["sha256"] / path.name
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        shutil.copy2(path, destination)
    return files.get(destination, asset["sha256"])


def snapshot_cache(cache, files, output, previous_payloads=None):
    """Publish immutable canonical transport copies for node-local cache files."""
    source = Path(cache["directory"])
    if not str(source).startswith("/mnt/data/"):
        return dict(
            cache,
            storage="shared",
            original_directory=str(source),
            canonical_directory=str(source),
        )
    identity = hashlib.sha256(
        json.dumps(cache["files"], sort_keys=True).encode()
    ).hexdigest()
    root = output / "cache-mirrors" / identity
    destination = root / "cache" / source.name
    destination.mkdir(parents=True, exist_ok=True)
    (root / "workflow").mkdir(exist_ok=True)
    (root / "blueprint").mkdir(exist_ok=True)
    shutil.copy2(source.parents[1] / "workflow/init.yaml", root / "workflow/init.yaml")
    blobs = output / "cache-payloads"
    blobs.mkdir(exist_ok=True)
    for name, digest in cache["files"].items():
        if name in {"vcfcache_annotated.bcf", "vcfcache_annotated.bcf.csi"}:
            blob = blobs / digest
            if not blob.exists():
                previous = previous_payloads / digest if previous_payloads else None
                if previous and previous.is_file() and sha256sum(previous) == digest:
                    # Same shared filesystem: preserve and reuse verified bytes
                    # rather than recopying GB payloads after a preparation stop.
                    blob.hardlink_to(previous)
                    print(f"Reusing verified canonical payload: {digest}", flush=True)
                else:
                    temporary = blobs / (digest + ".partial")
                    print(
                        f"Staging canonical cache payload: {source / name}", flush=True
                    )
                    shutil.copy2(source / name, temporary)
                    if sha256sum(temporary) != digest:
                        raise RuntimeError("canonical cache transport SHA mismatch")
                    temporary.rename(blob)
            if not (destination / name).exists():
                (destination / name).symlink_to(blob)
        else:
            shutil.copy2(source / name, destination / name)
        files.get(destination / name, digest)
    return dict(
        cache,
        storage="node_local",
        original_directory=str(source),
        canonical_directory=str(destination),
    )


def node_cache(cache, node_root=Path("/mnt/data/tmp/vcfcache-final-cache-replicas")):
    """Use a verified original replica, or install the canonical bytes locally."""
    if cache.get("storage") != "node_local":
        return Path(cache["directory"])
    identity = hashlib.sha256(
        json.dumps(cache["files"], sort_keys=True).encode()
    ).hexdigest()
    root = node_root / identity
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        ready = root / "VERIFIED.json"
        if ready.exists():
            previous = json.loads(ready.read_text())
            if previous["expected"] == cache["files"]:
                directory = Path(previous["directory"])
                if all(
                    (directory / name).is_file()
                    and (directory / name).stat().st_size == stat["bytes"]
                    and (directory / name).stat().st_mtime_ns == stat["mtime_ns"]
                    for name, stat in previous["metadata"].items()
                ):
                    return directory
        original = Path(cache["original_directory"])
        valid = all(
            (original / name).is_file() and sha256sum(original / name) == digest
            for name, digest in cache["files"].items()
        )
        if valid:
            destination = original
        else:
            canonical = Path(cache["canonical_directory"])
            installed = root / "installed"
            if installed.exists():
                raise RuntimeError(
                    "previously installed private cache changed; inspect before explicit replacement"
                )
            free = shutil.disk_usage(root).free
            needed = sum((canonical / name).stat().st_size for name in cache["files"])
            if free < needed + 25 * 1024**3:
                raise RuntimeError(
                    "insufficient node-local scratch for canonical cache plus runtime reserve"
                )
            with tempfile.TemporaryDirectory(prefix="install-", dir=root) as temporary:
                stage = Path(temporary)
                destination = stage / "cache" / original.name
                destination.mkdir(parents=True)
                (stage / "workflow").mkdir()
                (stage / "blueprint").mkdir()
                shutil.copy2(
                    canonical.parents[1] / "workflow/init.yaml",
                    stage / "workflow/init.yaml",
                )
                for name, digest in cache["files"].items():
                    shutil.copy2(canonical / name, destination / name)
                    if sha256sum(destination / name) != digest:
                        raise RuntimeError(
                            "node-local canonical cache installation SHA mismatch"
                        )
                stage.rename(installed)
            destination = installed / "cache" / original.name
        metadata = {
            name: dict(
                bytes=(destination / name).stat().st_size,
                mtime_ns=(destination / name).stat().st_mtime_ns,
            )
            for name in cache["files"]
        }
        write_json_atomic(
            ready,
            dict(
                directory=str(destination),
                expected=cache["files"],
                metadata=metadata,
                reused_original=valid,
            ),
        )
        return destination
