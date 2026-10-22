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
If bindgen cannot find `stddef.h` on a partially installed Clang toolchain,
install the matching Clang resource headers; as a local GCC fallback set
`BINDGEN_EXTRA_CLANG_ARGS="-I$(gcc -print-file-name=include)"` for the build.

## Parallelism and acceleration

`params.yaml: threads` is the shared configured CPU budget, further bounded by
process affinity. `--engine-workers 0` automatically chooses up to eight contig
workers, with two foreground pipeline processes budgeted per worker. Explicit
`--engine-workers N` overrides that choice within the same budget. Unused slots
are assigned to BGZF compression. With fewer than two configured CPUs a
single pipeline still has multiple cooperating processes; use taskset or a
scheduler cgroup when a strict process-wide quota is required.

Normalization/lookup/routing and final reconstruction are contig-parallel.
Adjacent contigs are grouped by indexed record count to reduce process starts;
both annotation lookups are restricted to the shard's indexed contig regions.
Annotation runs **once**, on all misses in input contig order: no repeated
annotator startup or changed recipe batch partitioning. Contigs are never split
across workers, avoiding indel-boundary duplication and cross-chunk matching.
A single-contig input therefore cannot exploit contig parallelism.

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
The manifest restores even the original order of records sharing a coordinate.
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
