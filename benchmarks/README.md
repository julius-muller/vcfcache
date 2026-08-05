# VCFcache publication benchmarks

This directory contains the reproducible preparation and execution plan for the
VCFcache publication benchmarks. Large datasets and results are kept outside Git
under `/mnt/data/vcfcache_benchmarks`.

The small, non-publication fastVEP diagnostic is documented separately in
[fastvep_pilot/README.md](fastvep_pilot/README.md). It runs on ITCCcloud_dev and
does not consume production Slurm workers.

Publication and repository plots are rendered exclusively in R with ggplot2.
See [figures/README.md](figures/README.md) for immutable source snapshots,
current drafts, and the refresh workflow. The older Python SVG functions are
retained only as diagnostic/compatibility utilities and are not publication
plotting backends.

## Documents

- [MATERIALS_AND_METHODS.md](MATERIALS_AND_METHODS.md): publication-ready
  methods draft, with explicit placeholders for the final full-cohort run.
- [SOURCE_PROVENANCE.md](SOURCE_PROVENANCE.md): source releases, citations,
  checksum chain, data-use notes, and submission-time provenance gate.
- [manifests/](manifests/): tracked machine-readable source, sample-selection,
  interval, and bundled-Zenodo-cache snapshots for peer review and deposition.
- [PLAN.md](PLAN.md): implementation phases and publication figure specification.
- [DATA_SETUP.md](DATA_SETUP.md): exact public-data layout, transformations, and QC.
- [INFRASTRUCTURE.md](INFRASTRUCTURE.md): measured VM capacity and scaling decision.
- [STORAGE_LAYOUT.md](STORAGE_LAYOUT.md): live dev/prod storage topology, data
  ownership, cache locations, result paths, and remaining filesystem capacity.
- [PILOT_RESULTS.md](PILOT_RESULTS.md): frozen paired-pilot results and scale-up gate.
- [ASSAY_DATA_RESULTS.md](ASSAY_DATA_RESULTS.md): completed HPRC/WES/panel counts,
  checksums, storage, and QC gate.
- [CACHE_STRATEGY_COMPARISON.md](CACHE_STRATEGY_COMPARISON.md): emergency
  validity assessment, runtime model, replacement assay/cache matrix, and the
  quarantined exploratory strategy comparison.
- [EXTERNAL_COHORT_SIMILARITY.md](EXTERNAL_COHORT_SIMILARITY.md): exact
  allele-level Jaccard and directional sharing checks for representative
  external WGS pairs.
- [figures/](figures/): provenance-tracked source snapshots and the canonical
  R/ggplot2 plotting workflow.

## Prepare the public VCF cohort

```bash
.venv/bin/python benchmarks/prepare_data.py preflight
.venv/bin/python benchmarks/prepare_data.py select
.venv/bin/python benchmarks/prepare_data.py smoke
.venv/bin/python benchmarks/prepare_data.py all
```

The default cohort contains 50 unrelated 1000 Genomes samples (10 from each
superpopulation, sex-balanced) and the seven GIAB v4.2.1 GRCh38 benchmark
genomes. Final inputs are sorted, BGZF-compressed, tabix-indexed, single-sample
VCFs covering autosomes 1–22.

All stages are resumable. Successfully verified downloads and completed
chromosome/sample outputs are skipped on subsequent invocations.

After QC, all files are also available through the flat, non-duplicating
directory `/mnt/data/vcfcache_benchmarks/samples/GRCh38/all`.

## Prepare the assay-size and HPRC R2 cohorts

The second preparation pipeline leaves the original 57-sample cohort immutable
and adds three assay/robustness cohorts:

- 20 population-stratified HPRC Release 2 graph-derived GRCh38 VCFs;
- 50 matched capture-like WES VCFs using the official Twist Human Core Exome
  hg38 targets padded by 125 bp, calibrated to approximately 80,000 variants;
- 50 matched strict-target WES controls using the unpadded Twist BED; and
- 50 matched small-panel VCFs covering the 84 ACMG SF v3.3 genes, defined as
  Ensembl 115 MANE Select coding sequence plus 20 bp on each side.

The word "independent" is intentionally not used here. The HPRC selection also
uses 1000 Genomes identities, while the WES and panel inputs are interval
subsets of the primary 1000 Genomes genomes. See the emergency validity
assessment in `CACHE_STRATEGY_COMPARISON.md` before interpreting bundled-cache
hit rates.

The matched WES and panel cohorts use chromosomes 1–22 and X. The official
1000 Genomes chrX callset is prepared as assay-only shards, so the established
autosomal WGS cohort does not change.

```bash
.venv/bin/python benchmarks/prepare_assay_data.py preflight
.venv/bin/python benchmarks/prepare_assay_data.py download --workers 3
.venv/bin/python benchmarks/prepare_assay_data.py select-hprc
.venv/bin/python benchmarks/prepare_assay_data.py regions
.venv/bin/python benchmarks/prepare_assay_data.py prepare-x
.venv/bin/python benchmarks/prepare_assay_data.py prepare-hprc
.venv/bin/python benchmarks/prepare_assay_data.py prepare-wes --workers 4
.venv/bin/python benchmarks/prepare_assay_data.py prepare-wes-targets --workers 4
.venv/bin/python benchmarks/prepare_assay_data.py prepare-panel --workers 4
.venv/bin/python benchmarks/prepare_assay_data.py qc
```

`prepare_assay_data.py all --workers 4` runs the same resumable sequence.
All final files remain ordinary sorted, BGZF-compressed, tabix-indexed,
single-sample `.vcf.gz` files; no BCF dataset is created.

The primary `wes_twist_core` cohort uses a transparent ±125-bp capture
footprint: 73,266–95,429 records per sample (median 77,056; mean 80,211). The
unpadded 33-Mb files are retained as `wes_twist_core_targets`, a mechanistic
target-only control rather than a claim about typical raw WES output.

## Tests

```bash
.venv/bin/python -m pytest tests/test_benchmark_prepare.py -q
.venv/bin/python -m pytest tests/test_benchmark_assay_prepare.py -q
.venv/bin/python -m pytest -q
VCFCACHE_BENCHMARK_NETWORK=1 \
  .venv/bin/python -m pytest tests/test_benchmark_prepare.py -m benchmark_network -q
```

## Paired publication pilot

The runner uses the local gnomAD GRCh38 cache selected by AF≥0.01 in at least
one `joint.freq` stratum and its immutable VEP 115.2 `--everything` recipe.
Cached and uncached commands are identical except for VCFcache's `--uncached`
switch.

```bash
.venv/bin/python benchmarks/run_pilot.py preflight
.venv/bin/python benchmarks/run_pilot.py all
```

Each run records the exact command, Git commit, environment, wall/CPU time,
maximum RSS, filesystem counters, VCFcache variant counts, and cache hit rate.
The final comparison streams both BCFs and requires identical variant keys,
input AF/AC/AN and GT values, and equivalent CSQ sets. It canonicalizes only
split-allele order within a locus and CSQ item order. Differences confined to
the VEP `HGNC_ID` CSQ field are counted and reported but do not fail the
comparison because VEP 115.2 assigns that field non-deterministically according
to buffering and batch composition (Ensembl VEP issue
[#1959](https://github.com/Ensembl/ensembl-vep/issues/1959)). Any difference in
another CSQ field still fails the comparison.

The single-sample pilot is feasibility evidence only. Before manuscript
submission, replace all square-bracketed fields in
`MATERIALS_AND_METHODS.md` from the full Slurm archive and complete the
submission-time gate in `SOURCE_PROVENANCE.md`.

## Bundled-cache strategy gate

Publication-facing public-cache comparisons accept only the two GRCh38 VEP
115.2 `--everything` cache aliases shipped by VCFcache and published as
annotated Zenodo bundles. The runner downloads and MD5-verifies them before use:

```bash
.venv/bin/python benchmarks/run_strategy_comparison.py fetch-bundled
```

An unpublished local annotation or a blueprint-only Zenodo record fails only
the **ready-made-cache** gate. Official GRCh38 blueprints are valid source
artifacts for locally building alternative annotation scenarios and have their
own checksum-frozen allow-list:

```bash
.venv/bin/python benchmarks/run_strategy_comparison.py fetch-blueprints
```

Those derived outputs must be labelled locally built, with blueprint DOI,
annotation/parameter hashes, cache identity, and build time reported.
Cohort-derived caches are likewise reported separately as custom caches.

## Slurm cohort

Prepare and submit one campaign from the controller. This creates a one-pair
HG02079 validation smoke job followed by one paired measurement for each of 50
samples. Condition order is balanced across samples. A legacy diagnostic
warm-up manifest is frozen for compatibility, but is not submitted by the
default chain; samples, not technical reruns, are the independent units.

```bash
.venv/bin/python benchmarks/run_cohort.py prepare \
  --campaign-id primary-wgs-<commit>
.venv/bin/python benchmarks/run_cohort.py submit-chain \
  --campaign-id primary-wgs-<commit> --concurrency 6
.venv/bin/python benchmarks/run_cohort.py status \
  --campaign-id primary-wgs-<commit>
.venv/bin/python benchmarks/run_cohort.py collect \
  --campaign-id primary-wgs-<commit>
.venv/bin/python benchmarks/collect_paired_benchmark.py \
  --campaign-root /mnt/data/slurm-results/campaigns/primary-wgs-<commit>
```

Each task runs cached and uncached modes in separate Slurm steps on worker-local
storage, requires semantic equivalence, and atomically archives the full pair
below `/results/campaigns/<campaign-id>/`. Manifests live on the same shared
export. The controller sees it at `/mnt/data/slurm-results`; workers mount it as
`/results`. Failed attempts are retained separately and can be resubmitted with
the `submit` subcommand and a sparse `--task-ids` array specification.

## Independent external-WGS cohort

The source-overlap upper-bound campaign is complemented by held-out public
genomes from KPGP, SGDP, and Harvard PGP. Source discovery, assembly-specific
header validation, deterministic allocation, normalization, QC, and provenance
are resumable:

```bash
.venv/bin/python benchmarks/prepare_external_wgs.py preflight
.venv/bin/python benchmarks/prepare_external_wgs.py catalog
.venv/bin/python benchmarks/prepare_external_wgs.py probe-pgp
.venv/bin/python benchmarks/prepare_external_wgs.py select
.venv/bin/python benchmarks/prepare_external_wgs.py download
.venv/bin/python benchmarks/prepare_external_wgs.py prepare
.venv/bin/python benchmarks/prepare_external_wgs.py qc
.venv/bin/python benchmarks/prepare_external_wgs.py screen-relatedness
for cohort in kpgp sgdp pgp; do
  .venv/bin/python benchmarks/prepare_external_wgs.py build-cache --cohort "$cohort"
done
```

If KING reports exactly one related pair, the deterministic recovery command
preserves the evaluation set when possible, substitutes the next eligible
sample from the same allocation stratum, archives the rejected screening work,
and records the decision in `manifests/relatedness_replacements.json`:

```bash
.venv/bin/python benchmarks/prepare_external_wgs.py replace-related
.venv/bin/python benchmarks/prepare_external_wgs.py download
.venv/bin/python benchmarks/prepare_external_wgs.py prepare
.venv/bin/python benchmarks/prepare_external_wgs.py qc
.venv/bin/python benchmarks/prepare_external_wgs.py screen-relatedness
```

A realistic relationship may instead be retained within the evaluation cohort.
`retain-related-evaluation` reverses the latest replacement and deterministically
promotes an unrelated evaluation genome to training, ensuring that a relative
never crosses the cohort-cache training/evaluation boundary.

KPGP and SGDP use open DDBJ/NIG GRCh38 autosomal gVCFs. PGP uses its larger
native GRCh37/hg19 single-sample subset; numeric contigs are renamed to UCSC
form before normalization, without coordinate liftover. Phenotype metadata is
never collected. The frozen design uses three cache-building genomes disjoint
from 20 KPGP, 20 SGDP, and 12 PGP evaluation genomes. PLINK2 KING screening is
performed separately per assembly. Relationships may be retained and reported
within the evaluation cohort, but kinship greater than 0.0884 across the
cache-training/evaluation boundary fails the design before any benchmark
campaign can be prepared.

Each external task runs one common uncached baseline plus the two verified,
assembly-matched bundled Zenodo caches (any-stratum AF >=10% and AF >=1%) and
the corresponding three-genome cohort cache. GRCh37 samples never use GRCh38
caches or vice versa. Condition order is balanced, every cached output is
semantically compared with the common baseline, and only the documented VEP
HGNC_ID discrepancy is ignored.

```bash
.venv/bin/python benchmarks/run_external_cohort.py prepare \
  --campaign-id external-wgs-<commit>
.venv/bin/python benchmarks/run_external_cohort.py submit-chain \
  --campaign-id external-wgs-<commit> --concurrency 6
.venv/bin/python benchmarks/run_external_cohort.py status \
  --campaign-id external-wgs-<commit>
.venv/bin/python benchmarks/run_external_cohort.py collect \
  --campaign-id external-wgs-<commit>
```

The ggplot2 figure workflow evaluates
`T_eff = T_cached + T_build/S` across cohort sizes 1 through 1000 and retains
the exact source grid beside the rendered output.

Campaigns prepared before the single-pass design amendment may already contain
a complete, validated first pass under the historical `warmup` phase. Collect
that phase directly instead of rerunning it:

```bash
.venv/bin/python benchmarks/run_external_cohort.py collect \
  --campaign-id external-wgs-<commit> --phase warmup
```

## Controlled pipeline-cost and hit-rate experiment

The mechanism panel uses one real capture-like WES input (HG02374, about 76,000
variants), four annotation costs, and deterministic self-caches targeting 50%,
80%, 90%, 95%, and 100% hits. The annotation costs are vanilla VEP, VEP
`--everything`, and vanilla VEP with a no-output plugin that pauses 5 or 20 ms
per transcript consequence. The plugin runs only on misses and cannot alter the
annotation schema. This is a deliberately compact 24-task experiment: one
uncached baseline per pipeline and one cached run per pipeline/hit-rate cell,
with no technical repeats.

Prepare the self-caches on the data-preparation VM, stage the resulting
`controlled_runtime` tree unchanged to every worker, then prepare and submit:

```bash
.venv/bin/python benchmarks/prepare_controlled_runtime.py
.venv/bin/python benchmarks/run_controlled_cohort.py prepare \
  --campaign-id controlled-runtime-<commit>
.venv/bin/python benchmarks/run_controlled_cohort.py submit \
  --campaign-id controlled-runtime-<commit> --concurrency 6
.venv/bin/python benchmarks/run_controlled_cohort.py collect \
  --campaign-id controlled-runtime-<commit>
```

The analysis fixes the miss-cost slope to the prespecified runtime equation and
fits only the nonnegative median lookup/preprocessing overhead. The final
controlled metrics and fitted values feed the R/ggplot2 mechanism panel and the
simple repository graphic so its overhead is empirical.

After the primary and assay campaigns are complete, collect their paired
metrics separately and combine them into the manuscript assay panel:

```bash
.venv/bin/python benchmarks/collect_paired_benchmark.py \
  --campaign-root /mnt/data/slurm-results/campaigns/assay-singlepass-<commit>
Rscript --vanilla benchmarks/figures/R/render_all.R \
  /path/to/frozen/source-snapshot /path/to/figure-output
```
