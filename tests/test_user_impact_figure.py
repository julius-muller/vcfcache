import pytest

from benchmarks.render_user_impact_figure import (
    DEFAULT_SCENARIOS,
    cached_seconds,
    human_time,
    impact_rows,
)


def test_cached_runtime_scaling_rule():
    assert cached_seconds(3600, 0.8, 0) == pytest.approx(720)
    assert cached_seconds(3600, 0.8, 30) == pytest.approx(750)
    assert cached_seconds(100, 0, 20) == 100


def test_impact_rows_cover_user_questions():
    rows = impact_rows(DEFAULT_SCENARIOS, (10, 60, 600), (2, 10, 100, 1000), 0, 0.8, 0)
    assert {row["panel"] for row in rows} == {
        "time_left",
        "pipeline_cost",
        "sample_scale",
    }
    scale = [row for row in rows if row["panel"] == "sample_scale"]
    assert [row["sample_count"] for row in scale] == [2, 10, 100, 1000]
    assert scale[-1]["saved_seconds"] == pytest.approx(100 * scale[1]["saved_seconds"])


def test_human_time_is_compact():
    assert human_time(30) == "30 sec"
    assert human_time(600) == "10 min"
    assert human_time(150) == "2.5 min"
    assert human_time(7200) == "2.0 h"
