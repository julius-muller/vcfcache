"""mismatch_report compares records in file order, also within a locus."""

import gzip
import json

from benchmarks.gather_engine import cohort


def _fake_loci(files):
    def loci(path):
        key, lines = None, []
        for line in files[path]:
            current = b"\t".join(line.split(b"\t")[:2])
            if current != key:
                if lines:
                    yield key, lines
                key, lines = current, []
            lines.append(line)
        if lines:
            yield key, lines

    return loci


def _rec(pos, alt, info):
    return f"chr1\t{pos}\t.\tA\t{alt}\t.\t.\t{info}\n".encode()


def _report(monkeypatch, tmp_path, candidate, reference):
    monkeypatch.setattr(cohort, "loci", _fake_loci({"c": candidate, "r": reference}))
    out = tmp_path / "m.jsonl.gz"
    result = cohort.mismatch_report("c", "r", out)
    with gzip.open(out, "rt") as stream:
        kept = [json.loads(line) for line in stream]
    return result, kept


def test_identical_files_report_nothing(monkeypatch, tmp_path):
    recs = [_rec(10, "C", "X=1;Y=2"), _rec(10, "G", "X=1"), _rec(20, "T", "X=3")]
    result, kept = _report(monkeypatch, tmp_path, recs, list(recs))
    assert result == dict(report_loci=0, report_loci_as_sets=0, report_loci_reordered=0)
    assert kept == []


def test_records_swapped_within_a_locus_are_reported(monkeypatch, tmp_path):
    a, b = _rec(10, "C", "X=1"), _rec(10, "G", "X=2")
    result, kept = _report(monkeypatch, tmp_path, [b, a], [a, b])
    assert result == dict(report_loci=1, report_loci_as_sets=0, report_loci_reordered=1)
    assert kept[0]["only_order_differs"] is True


def test_info_tag_order_alone_is_not_a_difference(monkeypatch, tmp_path):
    result, kept = _report(
        monkeypatch, tmp_path, [_rec(10, "C", "Y=2;X=1")], [_rec(10, "C", "X=1;Y=2")]
    )
    assert result["report_loci"] == 0
    assert result["report_loci_reordered"] == 0


def test_value_differences_are_reported(monkeypatch, tmp_path):
    result, kept = _report(
        monkeypatch, tmp_path, [_rec(10, "C", "X=1")], [_rec(10, "C", "X=2")]
    )
    assert result == dict(report_loci=1, report_loci_as_sets=1, report_loci_reordered=0)
    assert kept[0]["differs_as_sets"] is True
