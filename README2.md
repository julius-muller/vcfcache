[![DOI](https://zenodo.org/badge/947952659.svg)](https://zenodo.org/badge/latestdoi/947952659)
[![CI](https://github.com/julius-muller/vcfcache/actions/workflows/ci.yml/badge.svg)](https://github.com/julius-muller/vcfcache/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/julius-muller/vcfcache)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/vcfcache)](https://pypi.org/project/vcfcache/)
[![Cite](https://img.shields.io/badge/Cite-CITATION.cff-blue)](CITATION.cff)
[![codecov](https://codecov.io/github/julius-muller/vcfcache/graph/badge.svg?token=ELV3PZ6PNL)](https://codecov.io/github/julius-muller/vcfcache)

# VCFcache

VCFcache accelerates variant annotation by reusing annotations for variants that are already present in a local cache. Variants not found in the cache are annotated with your existing workflow and merged into a single output file.

VCFcache is designed to be a minimal wrapper around established annotation pipelines (e.g., VEP). It does not introduce a new annotation model; it reuses precomputed results where possible. It is genome‑agnostic and tool‑agnostic (VEP, SnpEff, ANNOVAR, custom scripts).

## When VCFcache helps

VCFcache is useful when you repeatedly annotate many samples with a stable pipeline (same reference build, same tool version, same options/plugins). Runtime improvements depend on the cache hit rate and I/O.

## Key properties

- **Drop-in integration:** your annotation command stays your annotation command; you place it into a simple `annotation.yaml`.
- **Cache reuse + fallback:** cache hits are looked up; cache misses are annotated using the same recipe.
- **BCF-native:** VCFcache reads and writes **BCF**. If your pipeline uses VCF/VCF.gz, convert at the edges with `bcftools view`.
- **Reproducibility boundaries:** cache correctness depends on using the same reference build and the same annotation recipe/tool versions used to create the cache.

---

## Quick start (recommended): use a published human cache (VEP)

VCFcache can discover and download published blueprints/caches from Zenodo.

### 1) Install

```bash
uv venv .venv && source .venv/bin/activate
uv pip install vcfcache
````

### 2) List available published caches

```bash
vcfcache list caches
vcfcache list blueprints
```

Tip: set a persistent download location (recommended for large caches):

```bash
export VCFCACHE_DIR=/path/to/large/disk
```

### 3) Run annotation using a published cache

After the cache is downloaded, annotate your input. VCFcache operates on BCF, so convert VCF/VCF.gz at the boundary.

**VCF.gz → VCFcache (BCF) → VCF.gz:**

```bash
bcftools view -Ob -o - sample.vcf.gz \
| vcfcache annotate -a /tmp/caches/<cache_name> -i - -o - \
| bcftools view -Oz -o sample.annotated.vcf.gz -W
```

Notes:

* Use `-i -` / `-o -` to stream from stdin / to stdout.
* If you already have BCF input, you can skip the first `bcftools view`.
* To keep BCF output, omit the final `bcftools view -Oz ...` step.

Optional: write runtime statistics (recommended for first evaluation):

```bash
vcfcache annotate -a /tmp/caches/<cache_name> -i input.bcf -o output.bcf --stats-dir stats/
```

---

## How to integrate VCFcache into an existing pipeline

Most users start with a working command such as:

```bash
vep ... --plugin X,my.bed --cache ... -i my.vcf -o my_annotated.vcf
```

With VCFcache you keep the toolchain, but you express the annotation step as a recipe that:

1. reads **BCF** from VCFcache,
2. converts to VCF if your tool requires it,
3. runs the annotator,
4. converts back to **BCF**.

### 1) Create an `annotation.yaml`

Minimal example for VEP (BCF in → VCF → VEP → BCF out):

```yaml
# annotation.yaml
annotation_cmd: |
  ${params.bcftools_cmd} view -Ov ${INPUT_BCF} | \
  ${params.annotation_tool_cmd} --offline --cache --vcf \
        --dir_cache ${params.vep_cache} \
        --assembly ${params.genome_build} \
        -i stdin -o stdout | \
  ${params.bcftools_cmd} view -Ob -o ${OUTPUT_BCF}
must_contain_info_tag: CSQ
required_tool_version: "115.2"
genome_build: "GRCh38"
```

Conventions:

* `INPUT_BCF` and `OUTPUT_BCF` are provided by VCFcache.
* Keep all tool options (plugins, custom files, flags) inside this command.
* If your command depends on local paths (FASTA, plugin data, cache dirs), put them in `params.yaml` and reference them as `${params.*}` here.

### 2) Provide environment-specific parameters (recommended)

```yaml
# params.yaml
genome_build: GRCh38
bcftools_cmd: bcftools
annotation_tool_cmd: vep
vep_cache: /path/to/vep/cache
threads: 8
temp_dir: /tmp
```

**Important:** any paths that differ between cache‑build time and annotation time must live in `params.yaml`.  
Example: instead of hardcoding `--plugin X,/path/to/data` in `annotation.yaml`, define:
```yaml
plugin_x_dir: /path/to/data
```
and reference it in `annotation.yaml` as:
```
--plugin X,${params.plugin_x_dir}
```

### 3) Run VCFcache using your recipe

```bash
vcfcache annotate -a /tmp/caches/<cache_name> \
  -y params.yaml \
  -i input.bcf -o output.bcf
```

Before running a cache you downloaded from someone else, use:
```bash
vcfcache annotate --requirements -a /tmp/caches/<cache_name>
```
This prints the expected tool versions, required `${params.*}` keys, and a fully substituted annotation command.

---

## Output equivalence and verification

VCFcache does not modify the annotation logic of your recipe:

* For **cache hits**, it reuses the stored annotations from the cache.
* For **cache misses**, it runs your `annotation.yaml` command.

To validate equivalence for your environment and options, run a small sample once with caching disabled (or using an empty cache), once with caching enabled, and compare outputs at the variant/INFO level using `bcftools` (and, if needed, normalized representations).

Practical recommendation:

* Start with a small VCF/BCF (e.g., a few thousand variants).
* Use a fixed reference build and deterministic tool settings.
* Compare post-processed outputs (normalized, sorted, same compression) rather than relying on byte-identical files.

---

## Requirements and constraints

* **Reference build must match** (e.g., GRCh37 vs GRCh38). Cache hits require matching contigs and coordinates.
* **Contig naming must match** (e.g., `chr1` vs `1`). Convert contigs at the boundary if needed.
* **BCF I/O:** VCFcache operates on BCF; use `bcftools view` to bridge tools that only read/write VCF.
* **Recipe compatibility:** cached annotations are only valid for the same annotation recipe/tool version used to create the cache.

---

## Building your own cache (optional)

If a published cache does not match your needs (different VEP flags, plugins, transcript set, or reference), you can:

1. initialize a blueprint (variant set),
2. build an annotated cache using your `annotation.yaml`.

The CLI shows the exact commands for blueprint download and cache build:

```bash
vcfcache list blueprints
vcfcache blueprint-init --doi <DOI> -o <output_dir>

vcfcache cache-build --doi <DOI> -a annotation.yaml -n <name>
```

---

## Links

* Documentation and concepts: `WIKI.md`
* Releases and published artifacts: Zenodo records shown by `vcfcache list …`
* Issues: GitHub issue tracker
