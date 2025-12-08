from pathlib import Path

from vcfcache.manifest import find_alias, load_manifest


def test_manifest_load_and_find_alias(tmp_path):
    manifest_file = tmp_path / "m.yaml"
    manifest_file.write_text(
        """
- alias: cache-hg38-gnomad-4.1wgs-AF0100-vep-115.2-basic
  type: cache
  version: "0.3.0"
  genome: hg38
  source: gnomad
  release: 4.1wgs
  filt: AF0100
  tool: vep
  tool_version: "115.2"
  preset: basic
  doi: 10.5281/zenodo.12345
  updated_at: 2025-01-01
  md5: abc
"""
    )

    entries = load_manifest(str(manifest_file))
    hit = find_alias(entries, "cache-hg38-gnomad-4.1wgs-AF0100-vep-115.2-basic")
    assert hit is not None
    assert hit.doi == "10.5281/zenodo.12345"
