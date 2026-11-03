"""Command-wall attribution is a critical-path approximation, never CPU summing."""

import pytest

from benchmarks.python_adaptive_results.log_headroom import components


def test_parallel_shards_use_maximum_not_sum():
    text = """Step 0/4:
Command completed in 1.0s: bcftools norm
Step 1/4:
Command completed in 2.0s: bcftools annotate
Command completed in 3.0s: bcftools annotate
Step 2/4:
Command completed in 4.0s: bcftools filter
Command completed in 5.0s: bcftools filter
Command completed in 6.0s: bcftools concat
Step 3/4:
Command completed in 7.0s: bcftools view
Step 4/4:
Command completed in 8.0s: bcftools annotate
Command completed in 9.0s: bcftools annotate
Command completed in 10.0s: bcftools concat
Command completed in 11.0s: bcftools view
"""
    result = components(text)
    assert result["cache_overlay_seconds"] == 3
    assert result["select_misses_seconds"] == 11
    assert result["merge_overlay_seconds"] == 9
    assert result["annotate_misses_seconds"] == 7
    assert result["final_write_index_seconds"] == 11


def test_incomplete_log_is_rejected():
    with pytest.raises(AssertionError, match="incomplete log"):
        components("Step 0/4:\nCommand completed in 1.0s: bcftools norm")
