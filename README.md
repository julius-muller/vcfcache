[![DOI](https://zenodo.org/badge/947952659.svg)](https://zenodo.org/badge/latestdoi/947952659)
[![CI](https://github.com/julius-muller/vcfcache/actions/workflows/ci.yml/badge.svg)](https://github.com/julius-muller/vcfcache/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/julius-muller/vcfcache)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/vcfcache)](https://pypi.org/project/vcfcache/)
[![Cite](https://img.shields.io/badge/Cite-CITATION.cff-blue)](CITATION.cff)
[![codecov](https://codecov.io/github/julius-muller/vcfcache/graph/badge.svg?token=ELV3PZ6PNL)](https://codecov.io/github/julius-muller/vcfcache)

# VCFcache – cache once, annotate fast

VCFcache accelerates variant annotation by caching annotations for a user-defined variant set and reusing them across runs, so only novel variants are processed at runtime.

### When VCFcache helps

VCFcache is useful when you either (a) repeatedly annotate many samples with a stable pipeline, or (b) want to quickly apply common annotations (e.g., VEP --everything) to a large VCF/BCF. Speedup depends on the cache hit rate of the input sample and the per-variant cost of the original annotation pipeline.

### Key properties

* **Drop-in integration:** keep your existing annotation command; place it into a simple `annotation.yaml` and run `vcfcache annotate`.
* **Cache reuse with automatic fallback:** cache hits are reused; cache misses are annotated with your configured command and merged into one output.
* **Genome- and tool-agnostic:** works with arbitrary reference builds and organisms, and with any annotator or pipeline that can be expressed as a command (stdin/stdout or file-based).
* **Pre-built caches available:** published caches built from public aggregation resources and annotation tools can be downloaded and used immediately (currently hg19/hg38 annotated with VEP --everything; additional configurations can be generated on request); the same tooling can generate highly efficient custom caches for your specific options and datasets.
* **BCF-native I/O:** VCFcache reads and writes **BCF** for performance and indexing; use `bcftools view` to convert VCF/VCF.gz at the boundaries.

**Important**: to use a prebuilt cache, you must have the same annotation tool (and compatible version) installed locally.

See [WIKI.md](WIKI.md) for full documentation, performance notes, and cache distribution via Zenodo.

---

## Quick Start

### Installation

```bash
# Via pip (requires Python 3.11+ and bcftools >= 1.20)
uv pip install vcfcache

# Via Docker (includes bcftools)
docker pull ghcr.io/julius-muller/vcfcache:latest

# Via Apptainer
apptainer exec docker://ghcr.io/julius-muller/vcfcache:latest vcfcache --help
```

See [WIKI.md - Section 2 (Quick Start)](WIKI.md#2-quick-start) for development installation and troubleshooting.

---

## Minimal example

This minimal example shows the complete workflow using a public cache.

### 1. List available caches

```bash
vcfcache list caches
```

### 2. Pull a specific cache

```bash
# Cache auto-downloads on first use, or download explicitly:
vcfcache cache-build --doi 10.5281/zenodo.18189447
```

### 3. Check cache requirements

```bash
vcfcache annotate --requirements -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic
```

This shows:
- Required annotation tool version (e.g., VEP 115.2)
- Required params (reference cache paths, etc.)
- The exact annotation command that will run

**Critical**: Install the exact tool version shown here to ensure compatibility with the cache annotations.

### 4. Annotate your sample

```bash
vcfcache annotate \
  -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
  --vcf sample.bcf \
  --output sample_annotated.bcf \
  --stats-dir ./results
```

**Input format**: The `--vcf` flag accepts VCF, VCF.gz, BCF, or stdin. VCFcache handles conversion automatically:
```bash
vcfcache annotate -a <cache-alias> --vcf sample.vcf.gz --output out.bcf     # VCF.gz input
vcfcache annotate -a <cache-alias> --vcf sample.vcf --output out.bcf        # VCF input
vcfcache annotate -a <cache-alias> --vcf sample.bcf --output out.bcf        # BCF input
bcftools view sample.vcf | vcfcache annotate -a <cache-alias> -i - -o -     # stdin/stdout
```

See [WIKI.md - Section 7 (Using a cache to annotate samples)](WIKI.md#7-using-a-cache-to-annotate-samples) for all annotation options.

---

## Building Your Own Cache

If you need different annotation settings (plugins, flags, tool version):

### From a public cache's blueprint
```bash
# Downloaded caches include the blueprint! Use it to build a cache variant:
# 1. Download cache (includes blueprint)
vcfcache cache-build --doi 10.5281/zenodo.18189447

# 2. Create annotation.yaml with your pipeline (see below)

# 3. Build cache from the downloaded blueprint
vcfcache cache-build \
  --db ~/.cache/vcfcache/caches/<cache-alias> \
  --name my_vep_cache \
  -a annotation.yaml \
  -y params.yaml
```

### From a public blueprint only
```bash
# 1. Download blueprint (smaller, no annotations)
vcfcache blueprint-init --doi <blueprint_DOI> -o ./cache_root

# 2. Create annotation.yaml with your pipeline (see below)

# 3. Build cache
vcfcache cache-build \
  --db ./cache_root \
  --name my_vep_cache \
  -a annotation.yaml \
  -y params.yaml
```

### From your own variants
```bash
# Use your cohort's common variants for maximum cache hit rate
vcfcache blueprint-init --vcf cohort_common.bcf --output ./cache_root -y params.yaml
vcfcache cache-build --db ./cache_root --name my_cache -a annotation.yaml -y params.yaml
```

See [WIKI.md - Section 6 (Building your own cache)](WIKI.md#6-building-your-own-cache-end-to-end) for the complete workflow including sharing via Zenodo.

---

## Setting Up annotation.yaml

The annotation.yaml defines your annotation pipeline **at cache-build time**. It is **immutable**: once baked into a cache, it cannot be changed without rebuilding.

### Start with your existing command

Begin with the exact command you already run:

```bash
vep --offline --cache --vcf \
  --dir_cache /data/vep/cache \
  --plugin ExACpLI,/data/plugins/ExACpLI.tsv.gz \
  --assembly GRCh38 \
  -i my.vcf.gz -o my_annotated.vcf.gz
```

### Adapt for VCFcache

Inside annotation.yaml, `${INPUT_BCF}` is always BCF (regardless of what format you provide at runtime). Since most annotation tools (like VEP) require VCF input, you need to convert it:

**Step 1**: Wrap for BCF↔VCF conversion:

```bash
# ▼ Convert BCF input to VCF for your tool
${params.bcftools_cmd} view -Ov ${INPUT_BCF} | \
  # ▼ Your existing VEP command (change files to stdin/stdout)
  ${params.annotation_tool_cmd} --offline --cache --vcf \  # ◄ Tool command (vep, docker run, etc.)
    --dir_cache /data/vep/cache \                          # ◄ Environment-specific path
    --plugin ExACpLI,/data/plugins/ExACpLI/ExACpLI.tsv.gz \ # ◄ Environment-specific path
    --assembly GRCh38 \
    -i stdin -o stdout | \                                 # ◄ Changed from -i/-o files
  # ▼ Convert VCF output back to BCF
  ${params.bcftools_cmd} view -Ob -o ${OUTPUT_BCF}
```

**What to change:**
- Replace `bcftools` → `${params.bcftools_cmd}`
- Replace `vep` → `${params.annotation_tool_cmd}` (allows docker, custom paths, etc.)
- Add `${params.bcftools_cmd} view -Ov ${INPUT_BCF} |` at the start (BCF → VCF)
- Add `| ${params.bcftools_cmd} view -Ob -o ${OUTPUT_BCF}` at the end (VCF → BCF)
- Change your tool's I/O from files to stdin/stdout
- Note which paths are environment-specific (marked with ◄)

**Step 2**: Final annotation.yaml with all substitutions:

```yaml
# annotation.yaml
annotation_cmd: |
  ${params.bcftools_cmd} view -Ov ${INPUT_BCF} | \
  ${params.annotation_tool_cmd} --offline --cache --vcf \
    --dir_cache ${params.vep_cache} \
    --plugin ExACpLI,${params.vep_plugins}/ExACpLI/ExACpLI.tsv.gz \
    --assembly GRCh38 \
    -i stdin -o stdout | \
  ${params.bcftools_cmd} view -Ob -o ${OUTPUT_BCF}

must_contain_info_tag: CSQ          # Required: VCFcache validates this tag exists
required_tool_version: "115.2"      # Required: VCFcache checks tool version match
genome_build: "GRCh38"              # Required: must match params.yaml
```

**Step 3**: Define environment-specific settings in params.yaml:

```yaml
# params.yaml
# Required fields:
genome_build: "GRCh38"              # Must match annotation.yaml
annotation_tool_cmd: "vep"          # Command to invoke annotation tool
bcftools_cmd: "bcftools"            # Command to invoke bcftools
temp_dir: "/tmp"                    # Directory for temporary files
threads: 8                          # Number of threads for bcftools

# Optional fields:
tool_version_command: "vep | awk '/ensembl-vep/ {print $NF}'"  # Command to check tool version (optional but recommended)

# Custom fields (for use in annotation.yaml):
vep_cache: "/data/vep/cache"        # Machine-specific path
vep_plugins: "/data/plugins"        # Machine-specific path
```

**About `tool_version_command`**: If provided, VCFcache runs this command and compares the output against `required_tool_version` in annotation.yaml. This validates that you have the correct tool version installed before annotating.

### Key rules

**Required fields in annotation.yaml**:
- `annotation_cmd`: Shell command to annotate variants
- `must_contain_info_tag`: INFO tag to validate (e.g., `CSQ` for VEP)
- `required_tool_version`: Tool version string (e.g., `"115.2"`)
- `genome_build`: Reference build (e.g., `GRCh38`) — **must match params.yaml**

**Required fields in params.yaml**:
- `annotation_tool_cmd`: Command to invoke annotation tool
- `bcftools_cmd`: Command to invoke bcftools
- `temp_dir`: Directory for temporary files
- `threads`: Number of threads for bcftools operations
- `genome_build`: Reference build (e.g., `GRCh38`) — **must match annotation.yaml**

**VCFcache-provided variables** (always use these):
- `${INPUT_BCF}`: Path to input BCF that VCFcache provides (the variants to annotate)
- `${OUTPUT_BCF}`: Path where annotated BCF should be written
- `${AUXILIARY_DIR}`: Directory for side outputs (HTML reports, logs, stats files, etc.)

**Environment-specific substitution** (when to use `${params.*}`):
- **Only needed if** you will use the cache on different machines with different paths
- **If substituting, substitute these**: tool paths, data directories, plugin files, reference paths
- **Don't substitute** annotation logic: flags like `--offline`, `--cache`, plugin names

**Immutability**:
- **At annotate time**, you can only override `params.yaml` (via `-y`)
- The `annotation.yaml` is frozen in the cache
- To change annotation logic, you must rebuild the cache

### Validation

Before annotating samples, check what the cache expects:

```bash
vcfcache annotate --requirements -a <cache-alias>
```

This shows:
- Required tool version and `${params.*}` keys
- The fully substituted annotation command
- Tool version validation results

See [WIKI.md - Section 8 (Configuration reference)](WIKI.md#8-configuration-reference-paramsyaml--annotationyaml) for full configuration details and advanced examples.

---

## Links

- **Full Documentation**: [WIKI.md](WIKI.md)
- **Performance Model**: [WIKI.md - Section 11](WIKI.md#11-performance-model-runtime-efficiency)
- **CLI Reference**: [WIKI.md - Section 10](WIKI.md#10-cli-reference-all-commands--flags)
- **Source**: https://github.com/julius-muller/vcfcache
- **Issues**: https://github.com/julius-muller/vcfcache/issues
- **Docker**: ghcr.io/julius-muller/vcfcache
