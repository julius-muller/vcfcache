#!/usr/bin/env python3
"""Correct pilot-only metadata in a freshly staged cache, never its source."""

import argparse
import json
import shutil
from pathlib import Path

import yaml

from benchmarks.adaptive_routing import digest


def main():
    """Move the redundant version-check key out of an annotation-only schema."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    args = parser.parse_args()
    recipe = args.cache / "annotation.yaml"
    config = yaml.safe_load(recipe.read_text())
    if "tool_version_command" not in config:
        return
    params = yaml.safe_load((args.cache / "params.snapshot.yaml").read_text())
    assert config["tool_version_command"] == params["tool_version_command"]
    backup = args.cache / "annotation.original.yaml"
    if backup.exists():
        raise FileExistsError(backup)
    shutil.copyfile(recipe, backup)
    del config["tool_version_command"]
    recipe.write_text(yaml.safe_dump(config, sort_keys=False))
    provenance = args.cache / "build_provenance.json"
    data = json.loads(provenance.read_text())
    data["original_recipe_sha256"] = data["recipe_sha256"]
    data["recipe_sha256"] = digest(recipe)
    data["metadata_repair"] = (
        "Removed redundant params-only tool_version_command; annotation_cmd and true build_seconds unchanged"
    )
    provenance.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
