"""Public cache manifest handling.

Defines a small schema for mapping cache aliases to Zenodo DOIs and metadata.
"""

import dataclasses
from pathlib import Path
from typing import List, Optional

import requests
import yaml


@dataclasses.dataclass
class CacheEntry:
    alias: str
    doi: str
    type: str  # blueprint | cache
    version: str
    genome: str
    source: str
    release: str
    filt: str
    tool: Optional[str] = None
    tool_version: Optional[str] = None
    preset: Optional[str] = None
    updated_at: str = ""
    md5: Optional[str] = None
    annotation_yaml_md5: Optional[str] = None


def load_manifest(path_or_url: str) -> List[CacheEntry]:
    if path_or_url.startswith("http://") or path_or_url.startswith("https://"):
        resp = requests.get(path_or_url, timeout=30)
        resp.raise_for_status()
        content = resp.text
    else:
        content = Path(path_or_url).read_text()

    data = yaml.safe_load(content) or []
    return [CacheEntry(**item) for item in data]


def find_alias(entries: List[CacheEntry], alias: str) -> Optional[CacheEntry]:
    return next((e for e in entries if e.alias == alias), None)


def format_manifest(entries: List[CacheEntry]) -> str:
    lines = ["alias\ttype\tversion\tgenome\tsource\trelease\tfilt\ttool\tdoi\tupdated"]
    for e in entries:
        lines.append(
            f"{e.alias}\t{e.type}\t{e.version}\t{e.genome}\t{e.source}\t{e.release}\t{e.filt}\t"
            f"{e.tool or ''}\t{e.doi}\t{e.updated_at}"
        )
    return "\n".join(lines)
