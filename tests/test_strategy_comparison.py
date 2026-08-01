from __future__ import annotations

from benchmarks.run_strategy_comparison import make_svg, select_samples


def test_strategy_selection_is_deterministic_and_disjoint(tmp_path):
    selection = tmp_path / "selected.tsv"
    selection.write_text(
        "sample\tpopulation\tsuperpopulation\tsex\tseed\n"
        + "".join(
            f"{pop}{index}\tP\t{pop}\tfemale\tseed\n"
            for pop in ("AFR", "AMR", "EAS", "EUR", "SAS")
            for index in range(3)
        )
    )
    first = select_samples(selection)
    second = select_samples(selection)
    assert first == second
    assert len(first["training"]) == 3
    assert len(first["evaluation"]) == 3
    assert {row["sample"] for row in first["training"]}.isdisjoint(
        row["sample"] for row in first["evaluation"]
    )
    assert [row["superpopulation"] for row in first["training"]] == first[
        "superpopulations"
    ]
    assert [row["superpopulation"] for row in first["evaluation"]] == first[
        "superpopulations"
    ]


def test_strategy_svg_contains_all_strategies_and_no_external_assets(tmp_path):
    rows = []
    strategies = (
        "gnomad_af_0.1",
        "gnomad_af_0.01",
        "gnomad_af_0.001",
        "cohort_3_genomes",
    )
    for strategy_index, strategy in enumerate(strategies):
        for sample_index in range(3):
            rows.append(
                {
                    "strategy": strategy,
                    "cache_hit_rate": 0.5 + 0.05 * strategy_index,
                    "speedup": 2 + strategy_index + sample_index / 10,
                }
            )
    output = tmp_path / "figure.svg"
    make_svg(output, rows)
    svg = output.read_text()
    assert svg.startswith('<svg xmlns="http://www.w3.org/2000/svg"')
    assert "gnomAD" in svg
    assert "3-genome" in svg
    assert "held-out genomes" in svg
    assert "href=" not in svg
