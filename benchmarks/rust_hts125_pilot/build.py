"""Build identical helper sources with bundled 1.19.1 and experimental HTSlib 1.25.

Only an isolated copy of hts-sys is replaced; Cargo registry/source worktrees
are never edited. Bindgen regenerates ABI bindings for each header version.
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tarfile
import tomllib
from pathlib import Path

ARCHIVE_SHA256 = "f96e583d2315998645bb56756f1c8f230c489512db2655b0f623198ebc8a2d99"


def digest(path):
    """Hash a build artifact without retaining its contents in memory."""
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    """Build and attest isolated ABI-matched helpers with pinned dependencies."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--crate", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=Path("native"))
    parser.add_argument("--cargo", default="cargo")
    parser.add_argument("--target-dir", type=Path)
    args = parser.parse_args()
    assert digest(args.archive) == ARCHIVE_SHA256
    args.root.mkdir(parents=True, exist_ok=False)
    crate = args.root / "hts-sys-125"
    shutil.copytree(args.crate, crate)
    (crate / "htslib").rename(crate / "original-htslib-119")
    # Retain the complete upstream release, including htscodecs, unmodified.
    with tarfile.open(args.archive) as archive:
        names = [
            m
            for m in archive.getmembers()
            if m.name.startswith("bcftools-1.25/htslib-1.25/")
        ]
        archive.extractall(args.root / "upstream", members=names, filter="data")
    shutil.copytree(args.root / "upstream/bcftools-1.25/htslib-1.25", crate / "htslib")
    # The hts-sys C source list/build options remain identical; fail if upstream
    # removed a source rather than silently building an incomplete library.
    environment = dict(os.environ)
    include = subprocess.check_output(
        ["gcc", "-print-file-name=include"], text=True
    ).strip()
    environment.setdefault("BINDGEN_EXTRA_CLANG_ARGS", "-I" + include)
    target = args.target_dir or args.root / "target"
    environment["CARGO_TARGET_DIR"] = str(target)
    builds = {}
    locks = {}
    for name, expected in (("old", "1.19.1"), ("latest", "1.25")):
        source = args.root / name
        source.mkdir()
        for file in ("Cargo.toml", "Cargo.lock"):
            shutil.copy2(args.source / file, source / file)
        shutil.copytree(args.source / "src", source / "src")
        command = [
            args.cargo,
            "build",
            "--offline",
            "--release",
            "-j4",
            "--manifest-path",
            str(source / "Cargo.toml"),
        ]
        if name == "old":
            command.append("--locked")
        else:
            command += ["--config", f'patch.crates-io.hts-sys.path="{crate}"']
        with (args.root / (name + "-build.log")).open("w") as log:
            subprocess.run(command, env=environment, stdout=log, stderr=log, check=True)
        binary = args.root / ("vcfcache-engine-" + name)
        shutil.copy2(target / "release/vcfcache-engine", binary)
        info = json.loads(
            subprocess.check_output([str(binary), "build-info"], text=True)
        )
        assert info["htslib_version"] == expected and info["libdeflate"]
        builds[name] = dict(
            binary=str(binary), sha256=digest(binary), info=info, command=command
        )
        locks[name] = tomllib.loads((source / "Cargo.lock").read_text())
        (args.root / "BUILD.json").write_text(
            json.dumps(
                dict(
                    builds=builds,
                    archive_sha256=ARCHIVE_SHA256,
                    source_hashes={
                        str(p.relative_to(args.source)): digest(p)
                        for p in sorted(args.source.rglob("*.rs"))
                    },
                    cargo_version=subprocess.check_output(
                        [args.cargo, "--version"], text=True
                    ).strip(),
                    dependency_lock_sha256=digest(args.source / "Cargo.lock"),
                ),
                indent=2,
            )
            + "\n"
        )
    # Dependency versions/features must not drift along with the C library.
    assert [p for p in locks["old"]["package"] if p["name"] != "hts-sys"] == [
        p for p in locks["latest"]["package"] if p["name"] != "hts-sys"
    ]
    assert (
        builds["old"]["info"]["htslib_features"]
        == builds["latest"]["info"]["htslib_features"]
    ), "Codec/build features differ"
    print(json.dumps(builds, indent=2))


if __name__ == "__main__":
    main()
