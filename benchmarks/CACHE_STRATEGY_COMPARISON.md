# Cache-strategy comparison

## Scope of the primary benchmark

The primary WGS Slurm campaign tests one annotation cache, not the entire
public-cache series:

- gnomAD v4.1 joint exome+genome sites, GRCh38;
- AF ≥ 0.01 (1%); 18,869,857 blueprint records;
- VEP 115.2, offline cache 115, `--everything`, `--hgvsg`, buffer 500,000,
  eight VEP forks and eight bcftools threads; and
- 3,950,678,136-byte annotated BCF, SHA-256
  `e9ae4a550c42a7cb1e68ccc1ef358c5e3558dab998c7bf833faa4f0571e1da7e`.

The filesystem suffix `af001` denotes the 0.01 threshold in these artifacts.
It must not be described as AF ≥ 0.001. The exact provenance is frozen in
[`SOURCE_PROVENANCE.md`](SOURCE_PROVENANCE.md).

## Small strategy experiment

The separate strategy experiment is descriptive and is not pooled with the
primary timing analysis. It compares identical VEP 115.2 recipes for:

1. gnomAD AF ≥ 0.1 (10%): 8,074,875 variants, public cache DOI
   [10.5281/zenodo.18189447](https://doi.org/10.5281/zenodo.18189447);
2. gnomAD AF ≥ 0.01 (1%): 18,869,857 variants, public cache DOI
   [10.5281/zenodo.18190046](https://doi.org/10.5281/zenodo.18190046);
3. gnomAD AF ≥ 0.001 (0.1%): 99,942,855 variants; the blueprint is public at
   [10.5281/zenodo.18190666](https://doi.org/10.5281/zenodo.18190666), while the
   locally retained 19.3-GB annotated cache is not currently listed as a
   production Zenodo cache; and
4. a custom cache made from the union of three primary-cohort genomes.

Training and evaluation identities are disjoint. The selection is frozen by
SHA-256 ranks using seed `vcfcache-paper-cache-strategy-v1`. The three selected
superpopulations and paired identities are:

| Superpopulation | Cache-training genome | Held-out evaluation genome |
|---|---|---|
| EAS | HG02050 | HG02028 |
| AMR | HG01148 | NA19759 |
| SAS | HG03850 | NA21124 |

This allocation is deterministic, prevents evaluating a sample against a
cache built from itself, and avoids choosing identities after observing
performance. One uncached run per evaluation genome supplies the common
baseline. Every cached output must pass the same semantic comparison against
that baseline, including the documented VEP issue #1959 treatment.

The figure has two panels: cache hit rate and end-to-end speedup versus
uncached annotation. Bars show the median and points show all three held-out
genomes. Cache record counts, cache bytes, raw wall times, hit rates, speedups,
and semantic comparison outcomes remain in the figure-source TSV.

## Reproduction

The preparation is resumable and deliberately runs outside the active Slurm
campaign on the data-preparation VM, where all cache artifacts already exist:

```bash
.venv/bin/python benchmarks/run_strategy_comparison.py prepare
.venv/bin/python benchmarks/run_strategy_comparison.py execute
.venv/bin/python benchmarks/run_strategy_comparison.py collect
```

`all` runs the same three stages in order. Outputs are written below
`/mnt/data/vcfcache_benchmarks/strategy_comparison`:

- `design.json`: frozen training/evaluation allocation;
- `custom_cache_build.json`: custom-cache build timing and completion state;
- `runs/`: raw commands, logs, metrics, BCFs, and semantic comparisons;
- `figure/cache_strategy.tsv`: publication figure source; and
- `figure/cache_strategy.svg`: editable vector figure.
