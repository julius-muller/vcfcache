from __future__ import annotations

from benchmarks.prepare_wgs_pipeline_spectrum import (
    DELAYS_US,
    delayed_annotation,
)

ANNOTATION = """genome_build: GRCh38
annotation_cmd: |
  vep \\
    --everything \\
    --stats_file ${AUXILIARY_DIR}/stats.html \\
    -o ${OUTPUT_BCF}
must_contain_info_tag: CSQ
"""


def test_wgs_spectrum_preserves_everything_and_adds_only_delay_plugin(tmp_path):
    rendered = delayed_annotation(ANNOTATION, 2_000, tmp_path / "plugins")
    assert "--everything" in rendered
    assert f"--dir_plugins {tmp_path / 'plugins'}" in rendered
    assert "--plugin SyntheticDelay,delay_us=2000" in rendered
    assert rendered.index("--plugin") < rendered.index("--stats_file")


def test_wgs_spectrum_has_six_monotonic_subday_virtual_loads():
    assert DELAYS_US == (500, 1_000, 2_000, 4_000, 7_000, 10_000)
