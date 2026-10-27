"""Untimed legacy-overlay audit reusing retained annotations, not reannotating."""

import argparse
import json
import shlex
import subprocess
import tempfile
from pathlib import Path

from benchmarks.engine_profile.run import read
from benchmarks.performance_sprint.manuscript import annotation_digest
from benchmarks.run_pilot import strict_semantic_compare, write_json_atomic


def main():
    """Reconstruct legacy overlay semantics without invoking an annotator again."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = read(args.campaign / "manifest.json")
    for identifier in ["d3a4e19d85523300", "3693ce45be7d4941"]:
        task = next(t for t in manifest["tasks"] if t["id"] == identifier)
        evidence = args.campaign / "attempts/a001" / identifier
        result_dir = args.output / identifier
        result_dir.mkdir()
        with tempfile.TemporaryDirectory(
            prefix="overlay-diagnostic-", dir="/mnt/data/tmp"
        ) as directory:
            work = Path(directory)
            source = work / "input.bcf"
            recipe = __import__("yaml").safe_load(
                (Path(task["cache"]["directory"]) / "annotation.yaml").read_text()
            )
            tag = recipe["must_contain_info_tag"]
            pipeline = [
                [
                    "bcftools",
                    "view",
                    task["input"]["path"],
                    "-R",
                    task["preprocess"]["bed"]["path"],
                    "-Ou",
                ],
                ["bcftools", "norm", "-m-", "-Ou"],
                ["bcftools", "view", "-e", 'ALT="*"', "-Ou"],
                [
                    "bcftools",
                    "norm",
                    "--rm-dup",
                    "exact",
                    "-Ob",
                    "-W",
                    "-o",
                    str(source),
                ],
            ]
            subprocess.run(
                [
                    "bash",
                    "-o",
                    "pipefail",
                    "-c",
                    " | ".join(shlex.join(p) for p in pipeline),
                ],
                check=True,
            )
            cache = str(Path(task["cache"]["directory"]) / "vcfcache_annotated.bcf")
            looked, misses, provided, merged = [
                work / f"{name}.bcf"
                for name in ["looked", "misses", "provided", "merged"]
            ]
            commands = [
                [
                    "bcftools",
                    "annotate",
                    "-a",
                    cache,
                    str(source),
                    "-c",
                    "INFO",
                    "-Ob",
                    "-W",
                    "--threads",
                    "8",
                    "-o",
                    str(looked),
                ],
                [
                    "bcftools",
                    "view",
                    "-i",
                    f'INFO/{tag}==""',
                    str(looked),
                    "-Ob",
                    "-W",
                    "-o",
                    str(misses),
                ],
                [
                    "bcftools",
                    "annotate",
                    "-a",
                    str(evidence / "output.bcf"),
                    str(misses),
                    "-c",
                    "INFO",
                    "-Ob",
                    "-W",
                    "-o",
                    str(provided),
                ],
                [
                    "bcftools",
                    "annotate",
                    "-a",
                    str(provided),
                    str(looked),
                    "-c",
                    "INFO",
                    "-Ob",
                    "-W",
                    "-o",
                    str(merged),
                ],
                [
                    "bcftools",
                    "view",
                    "-i",
                    f'INFO/{tag}!=""',
                    str(merged),
                    "-Ob",
                    "-W",
                    "-o",
                    str(result_dir / "legacy-composed.bcf"),
                ],
            ]
            for command in commands:
                subprocess.run(command, check=True)
            comparison = strict_semantic_compare(
                result_dir / "legacy-composed.bcf", evidence / "output.bcf"
            )
            comparison["legacy_misses"] = int(
                subprocess.check_output(
                    ["bcftools", "index", "-n", str(misses)], text=True
                )
            )
            comparison["new_counts"] = read(evidence / "statistics.json")[
                "variant_counts"
            ]
            comparison["legacy_digest"] = annotation_digest(
                result_dir / "legacy-composed.bcf", task["baseline"]
            )
            comparison["historical_digest"] = task["baseline"]["cached_digest"]
            write_json_atomic(result_dir / "comparison.json", comparison)
            print(
                json.dumps(
                    dict(
                        cell=identifier,
                        semantic_pass=comparison["semantic_pass"],
                        legacy_misses=comparison["legacy_misses"],
                        legacy_digest=comparison["legacy_digest"],
                        historical_digest=comparison["historical_digest"],
                        examples=comparison["examples"][:3],
                    )
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
