# Experimental fused engine

The Python CLI remains the public interface. `vcfcache-engine` is an optional
Rust/HTSlib subprocess used by `vcfcache annotate --engine fused`. Legacy remains
the default. Existing Zenodo caches and immutable recipes are unchanged.

## Build and use

Install a stable Rust toolchain, a C compiler, pkg-config, zlib development
headers and libclang development headers. Debian/Ubuntu build dependencies:

```sh
sudo apt-get install build-essential pkg-config zlib1g-dev libclang-dev
cargo build --release --locked --manifest-path native/Cargo.toml
export VCFCACHE_ENGINE="$(pwd)/native/target/release/vcfcache-engine"
vcfcache annotate -a /path/to/cache -i sample.bcf -o result.bcf --engine fused
```

Alternatively `cargo install --locked --path native` installs the executable on
PATH. No runtime compilation or Rust dependency is added to normal Python
installation. A repository checkout also discovers its release binary directly.
Adaptive routing requires helper version 0.2.0 or later. Rebuild after switching
branches; an explicitly installed older helper on PATH must be updated or
overridden with `VCFCACHE_ENGINE`. Frozen pilot snapshots retain their original
pre-release version labels and executable hashes.
If bindgen cannot find `stddef.h` on a partially installed Clang toolchain,
install the matching Clang resource headers; as a local GCC fallback set
`BINDGEN_EXTRA_CLANG_ARGS="-I$(gcc -print-file-name=include)"` for the build.

## Parallelism and acceleration

`params.yaml: threads` is the shared configured CPU budget, further bounded by
process affinity. `--engine-workers 0` automatically chooses up to eight contig
workers, initially budgeting four CPU slots per worker so threaded BGZF input
decoding is not starved by excessive process concurrency. Explicit
`--engine-workers N` overrides that choice within the same budget. Unused slots
are assigned to BGZF compression. During the fused workflow, the calling thread
and its newly spawned workers/descendants are restricted to the configured
number of already-allowed CPUs; the caller's affinity is restored even on
failure. Other pre-existing application threads are untouched. Affinity is a
core-concurrency limit, not a CPU-time quota; use a scheduler cgroup for quotas
and for restrictions covering CLI construction or unrelated application work.

After one-pass preparation, large inputs with estimated candidate cache counts
at most 1.25 times the input count instead use two slots per worker (up to eight
workers). This preserves foreground concurrency for near one-to-one dense
overlays; much larger caches retain decoder-priority allocation. CSI counts and
candidate-byte fractions determine this reproducible policy; it is not runtime
self-learning or a claim that all workloads scale linearly.
Inputs of at most 2,000 normalized records use one batched routing task by
default, avoiding per-shard process/filesystem overhead; available codec threads
can still accelerate its cache reads. Explicit worker overrides are retained.

Normalization is parallelized by input counts; cache lookup is scheduled by
conservative CSI byte estimates plus input work, not by sample counts alone.
During normalization a Rust helper writes indexed level-zero BCF and collects
16-KiB sample windows without a separate full-record planning scan.
`--engine-lookup auto` selects sparse windows below half the contig's estimated
byte footprint, otherwise contiguous access. The `regions` and `contigs`
overrides support diagnostics. Estimates union overlapping CSI chunks but
conservatively include the final BGZF block; they are **not observed I/O**.
Empty/missing-cache-contig inputs do not trigger a genome-wide pass. One worker
retains explicit selection.
Annotation runs **once**, on all misses in input contig order: no repeated
annotator startup or changed recipe batch partitioning. Large contigs may be
split into coordinate-disjoint tasks with POS-only selection; all records sharing
a position stay in the same task. Work is dynamically queued and outputs are
gathered in header/coordinate order, not completion order.

The native HTSlib build enables libdeflate's CPU-dispatched accelerated codecs;
the Python orchestrator does not parse or copy per-record data. Portable release
builds do not assume AVX2/AVX-512. SIMD benefits apply where supported by the
codec and CPU; there is no claim of a custom SIMD BCF matching kernel. HTSlib
and bcftools remain responsible for the binary format. Optional local-only
`RUSTFLAGS="-C target-cpu=native"` builds must not be distributed to other CPUs.

## Data path and compatibility

1. Normalize each input contig to minimally compressed indexed BCF. bcftools
   requires an indexed target when looking up annotations from BCF.
2. Stream the cache overlay directly into Rust. Write uncompressed hit/miss
   spools and a one-byte-per-record routing manifest; count records in-stream.
3. Collect/index the misses and invoke the original recipe once.
4. Overlay only the miss annotations onto their original sample records. Rust
   reconstructs input order, applies the final annotation filter and writes
   compressed shards directly using HTSlib/libdeflate and optional BGZF workers.
5. Header-checked naive concatenation copies compressed blocks without another
   whole-genome recompression. Build a final CSI, then publish the files.

This trades additional temporary disk space for fewer full-file scans and less
compression. Memory is bounded by workers and records at a shared coordinate,
not by total WGS size. Input and output sample metadata, GT/FORMAT, coordinates,
IDs, QUAL and FILTER are preserved according to legacy INFO-overlay semantics.
The manifest restores the cache-overlay stream's order within a coordinate;
bcftools itself may reorder distinct alleles at a shared position. Complete
record content and coordinate sorting, not byte-identical within-locus order
relative to another bcftools execution, are the compatibility contract.
This engine does not change the cacheability of sample-dependent or
context-dependent annotations; existing recipe requirements still apply.

Unsplit multiallelic/non-variant records, non-String required tags and duplicate
variant keys routed to both partitions cause an explicit legacy replay before
output is published. The effective engine and reason are recorded. Parsing,
subprocess and recipe failures are errors, **not** fallback triggers.
INFO overlays reject conflicting Type/Number declarations rather than coercing
values silently. This stricter invalid-input check is intentional.
Annotator output is indexed and sorted only if not indexable. Temporary spools survive
failures for diagnosis unless statistics are explicitly disabled; normal
successful workflow cleanup is unchanged.

`compare_stats.yaml` retains existing fields and adds requested/effective engine,
native version, worker/CPU settings and stage timings. Light-mode counters do
not rescan final records. Full audits remain optional. Diagnostics use stderr
so binary stdout stays clean. File output and its CSI are completed in a
destination-local staging directory before rename; the two renames are not a
filesystem transaction, so consumers must wait for successful completion.

## Derived-store experiment (not a production backend)

`vcfcache-engine store-build --cache CACHE.bcf --output NEW_DIRECTORY` explicitly
builds a read-only mmap prototype containing plain BCF payloads, full variant
keys and a locus index. It never overwrites an existing directory. The manifest
records the source checksum, conversion duration and storage amplification.
`store-query --store DIRECTORY --input NORMALIZED.bcf` emits raw candidate BCF;
pipe this through `bcftools view -Ob0 -W` before indexed use. All cache alleles at
a requested locus are retained so bcftools, rather than a new exact-key matcher,
continues to own padded-REF, allele-pair and INFO-overlay semantics. Hash matches
are verified against stored locus keys; collisions cannot identify variants.

These experimental commands are not automatically used by `annotate`. A store
must beat the **adaptive BCF path**, not the old full-cache scan, by at least 2x
in cache handling and 1.5x end-to-end on real fast-annotator WGS before a public
store backend is considered. Benchmark candidate construction, normalization and
indexing are all timed; one-time conversion is reported separately. Store files
must remain immutable while mapped. This prototype is an engineering test, not
a published cache format, DOI asset or new cache provenance category.

## Verification

```sh
cargo test --locked --manifest-path native/Cargo.toml
cargo clippy --locked --manifest-path native/Cargo.toml -- -D warnings
python -m pytest tests/test_fused_engine.py -q
docker build --target native-test -f docker/Dockerfile.vcfcache -t vcfcache:native-test .
docker run --rm vcfcache:native-test
```

Tests compare complete records and relevant headers against the legacy engine,
with no annotation-field exceptions. Performance results from the synthetic
runner are engineering diagnostics, not biological hit-rate estimates:

```sh
python -m benchmarks.fused_engine --output /scratch/fused-pilot --records 1000000
```

Scaling is constrained by contig sizes, storage/memory bandwidth, the remaining
lookup cost and the single recipe invocation. Linear scaling is a target, not
an assertion. No production default changes until correctness and measured
end-to-end benefit justify them.
