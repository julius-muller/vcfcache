# Cache-strategy comparison

> **Emergency review, 2026-08-01:** The original three-training-genome,
> mixed-superpopulation experiment below is retained as exploratory evidence
> only. It was stopped before the third evaluation genome completed and its
> output was marked `quarantined_exploratory`. It must not supply the
> publication figure. Slurm measured job 30 is held. The replacement design is
> driven by observed cache hit rate and includes more than ten panel, WES, and
> WGS inputs per strategy; see "Replacement assay and cache design".

## Scope of the primary benchmark

The primary WGS Slurm campaign tests one annotation cache, not the entire
public-cache series:

- gnomAD v4.1 joint exome+genome sites, GRCh38;
- AF ≥ 0.01 (1%) in at least one `joint.freq` stratum; 18,869,857 blueprint
  records;
- VEP 115.2, offline cache 115, `--everything`, `--hgvsg`, buffer 500,000,
  eight VEP forks and eight bcftools threads; and
- 3,950,678,136-byte annotated BCF, SHA-256
  `e9ae4a550c42a7cb1e68ccc1ef358c5e3558dab998c7bf833faa4f0571e1da7e`.

The filesystem suffix `af001` denotes the 0.01 threshold in these artifacts.
It must not be described as AF ≥ 0.001. The exact provenance is frozen in
[`SOURCE_PROVENANCE.md`](SOURCE_PROVENANCE.md).

This cache is broader than a cache selected by combined-population AF ≥ 0.01.
The exporter selected a site when *any* ancestry, sex, or subset frequency met
the threshold, but wrote combined `joint.freq[0]` AF into INFO. Publication
text and figure labels must therefore say **any-stratum AF ≥ 1%**, not merely
"gnomAD 1%". A separately rebuilt combined-population cache is not bundled and
is therefore ineligible for the bundled-cache benchmark.

## Quarantined small strategy experiment

The separate strategy experiment is descriptive and is not pooled with the
primary timing analysis. It compares identical VEP 115.2 recipes for:

1. gnomAD any-stratum AF ≥ 0.1 (10%): 8,074,875 variants, public cache DOI
   [10.5281/zenodo.18189447](https://doi.org/10.5281/zenodo.18189447);
2. gnomAD any-stratum AF ≥ 0.01 (1%): 18,869,857 variants, public cache DOI
   [10.5281/zenodo.18190046](https://doi.org/10.5281/zenodo.18190046);
3. gnomAD any-stratum AF ≥ 0.001 (0.1%): 99,942,855 variants; the blueprint is public at
   [10.5281/zenodo.18190666](https://doi.org/10.5281/zenodo.18190666), while the
   locally retained 19.3-GB annotated cache is not currently listed as a
   production Zenodo cache; and
4. a custom cache made from the union of three primary-cohort genomes.

The AF ≥ 0.1% entry is excluded from every replacement/publication analysis:
Zenodo record 18190666 is a blueprint, not a bundled annotated cache. Its local
annotation is not a downloadable VCFcache bundle and therefore does not meet
the cache-provenance gate.

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

The quarantined run was produced by commit `3039e57`; it must not be resumed.
The replacement runner uses a new output root, downloads the two allowed
bundled caches from production Zenodo, verifies the archive MD5 values, and
writes a provenance marker before preparation can start:

```bash
.venv/bin/python benchmarks/run_strategy_comparison.py fetch-bundled
.venv/bin/python benchmarks/run_strategy_comparison.py prepare
.venv/bin/python benchmarks/run_strategy_comparison.py execute
.venv/bin/python benchmarks/run_strategy_comparison.py collect
```

`all` runs download, preparation, execution, and collection in order. Public
strategies are restricted to these shipped aliases:

| Bundled alias | Zenodo DOI | Archive MD5 |
|---|---|---|
| `cache-gnomad-v4.1-GRCh38-joint-af01-vep115.2-e` | `10.5281/zenodo.18189447` | `088cf426461a51b77bdfcd5dcd2233f4` |
| `cache-gnomad-v4.1-GRCh38-joint-af001-vep115.2-e` | `10.5281/zenodo.18190046` | `3ac438461eac0cf42c75717156d7b2d4` |

The cohort-derived cache remains a distinct custom-cache strategy and is never
labelled bundled. Replacement outputs are written below
`/mnt/data/vcfcache_benchmarks/strategy_comparison_zenodo_v1`. Bundles and
retained archives live below
`/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches`.

Quarantined exploratory outputs remain below
`/mnt/data/vcfcache_benchmarks/strategy_comparison`:

- `design.json`: frozen training/evaluation allocation;
- `custom_cache_build.json`: custom-cache build timing and completion state;
- `runs/`: raw commands, logs, metrics, BCFs, and semantic comparisons;
- `figure/cache_strategy.tsv`: publication figure source; and
- `figure/cache_strategy.svg`: editable vector figure.

## Emergency validity assessment

The 50 read-called WGS inputs are NYGC 30x 1000 Genomes samples. The 20 HPRC
robustness inputs were also deliberately selected by intersecting HPRC with
1000 Genomes identities. gnomAD v3.1 added approximately 2,500 1000 Genomes
samples, and those genomes carry into gnomAD v4.1. Therefore neither local WGS
set is independent of the public-cache source universe. The existing gnomAD
results can calibrate lookup overhead and the relationship between miss count
and runtime, but cannot estimate the absolute cache hit rate expected for an
independent clinical WGS cohort.

The two completed exploratory genomes confirm that the hit-rate denominator is
not the problem: it is the number of output variants minus the variants sent to
VEP, divided by total output variants. The completed points fit

`T_cached = 7.62 min + 125.82 min * (1 - f)`

with R² approximately 0.996 and median uncached time approximately 127 min.
This is useful calibration evidence for the runtime model, despite the
non-independent samples.

## Replacement assay and cache design

The publication benchmark separates two questions:

1. **Runtime mechanics.** Fit lookup/preprocessing overhead and uncached
   per-variant annotation cost from real runs. Plot relative runtime against
   the *observed* hit rate rather than using training-cohort size as a proxy.
2. **Expected hit rate.** Measure the hit-rate distribution in an appropriate
   cohort, separately for public and cohort-derived caches. Do not infer this
   distribution from the overlapping 1000 Genomes/gnomAD comparison.

For each assay, run at least 12 held-out inputs through the same VEP 115.2
`--everything` recipe. Record total variants, hits, misses, lookup overhead,
cached wall time, uncached wall time, semantic equivalence, and cache-build
time. Report cached relative runtime as

`R_cached = T_cached / T_uncached ≈ T_overhead / T_uncached + (1 - f)`

and speedup as `1 / R_cached`. Report build-amortized relative runtime for
several use counts `S` as

`R_eff(S) = (T_cached + T_build / S) / T_uncached`.

The empirical matrix is:

| Assay | Available now on every worker | Permitted use |
|---|---:|---|
| ACMG SF panel | 50 | Assay scaling and cohort-cache evaluation |
| Capture-like WES | 50 | Assay scaling and cohort-cache evaluation; in-silico capture label required |
| NYGC 30x WGS | 50 | Runtime calibration and cohort-cache evaluation only |
| HPRC R2 WGS | 20 | Separate graph/assembly robustness analysis only |

For cohort-derived caches, use deterministic disjoint training and held-out
sets and report the actual hit-rate distribution. Exact training sizes such as
5, 10, or 100 are secondary metadata, not the x-axis; select enough frozen
training states to span the observed WGS range and show how build cost
amortizes.

For the bundled-cache WGS arm, acquire at least 12 read-called genomes whose
identities are not present in gnomAD. Until that provenance gate passes, the
1000 Genomes public-cache points must be labelled **source-overlap upper
bound** and excluded from the headline WGS hit-rate estimate. The matched WES
and panel inputs remain valuable for controlled assay-size comparisons, but
because they derive from the same genomes they do not repair the WGS source
overlap.
