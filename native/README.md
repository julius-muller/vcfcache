# Experimental native engines

The Python CLI remains the public interface. `vcfcache-engine` is an optional
Rust/HTSlib subprocess used by `vcfcache annotate --engine fused`. Legacy remains
the default. Existing Zenodo caches and immutable recipes are unchanged.

## Resident WES prototype (`--engine native`)

Requires helper **0.3.0**. This opt-in engine reads standard BCF/CSI directly;
there is no new storage format, service or annotator-specific integration.

```sh
vcfcache annotate -a /path/to/cache -i sample.bcf -o result.bcf \
  --engine native --engine-memory-mb 512
```

bcftools retains multiallelic splitting. Its normalized stream enters one
resident helper through a FIFO rather than an indexed sample spool. Independent
cache readers batch occupied 16-KiB windows, skip unpacking unrelated candidates,
and share a bounded BGZF decoder pool. Default lookup workers are
`min(8, max(1, available_cpus/2))`; remaining slots decode. Explicit worker counts
cannot exceed the CPU budget. Output compression uses the released budget.
`auto`/`regions` select occupied windows; `contigs` is a diagnostic whole-contig
override. These native policies differ from the older fused planner.

Original records and cached INFO payloads stay in memory. One indexed miss BCF
is passed to the **unchanged external recipe once**. Recipe annotations are
merged while consuming original records into one compressed, indexed output;
there are no normal-path hit spools, shard assembly or final recompression.
INFO is overlaid; original sample FORMAT/GT and other record fields remain.

The memory option bounds conservatively estimated retained-record buffers, not
whole-process RSS or the external annotator. Reader/header/codec allocations and
one recipe-output locus are additional. An indexed record-count lower bound can
reject WGS immediately; actual buffered growth is also checked before annotation.
Overflow, ambiguous duplicate keys, unsplit multiallelic records and potentially
compatible padded REF matches replay fused before invoking the recipe. No input
or annotation format is silently changed to force a faster path.

If recipe output needs non-exact pairing, duplicate resolution or contig-order
translation, a compatibility INFO overlay reuses its completed output; the
recipe is not rerun. Parse/subprocess errors are failures, not silent fallbacks.
Output is staged before publication; its CSI is checked against the first
written locus. Diagnostics, buffer estimates, candidates, queries and phase
timings are recorded in statistics. `merge-and-write` includes the streaming
recipe overlay and output encoding, not just compression.

See [pilot results and limitations](../benchmarks/native_bcf/RESULTS.md).

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
