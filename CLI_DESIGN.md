# CLI evolution: wrap an existing annotation command

Status: design proposal, not an implemented interface. Recorded 2026-10-22.
The existing recipe-based CLI remains supported. This document captures interface
ideas for a future feature branch based on the fused-routing branch; it does not
change runtime behaviour or promise support for the example commands below.

## Goal

Make caching as easy to adopt as prefixing an existing annotation command:

```bash
# Existing workflow
myannotator input.vcf -o output.bcf -plugin XX ...

# Desired workflow, for a supported annotator
vcfcache myannotator input.vcf -o output.bcf -plugin XX ...
```

The wrapper should discover an exactly compatible cache, reuse its annotations,
annotate misses and assemble the requested output. An agent may help configure a
pipeline, but should not be necessary for routine, safe use.

## Two entry points, one execution backend

Offer a convenient shorthand for supported annotators and an unambiguous form
that separates wrapper options from the annotator command:

```bash
# Proposed shorthand
vcfcache vep -i input.vcf -o output.vcf --plugin XX ...

# Proposed explicit form: wrapper options before --, annotator arguments after it
vcfcache run --threads 16 -- vep -i input.vcf -o output.vcf --plugin XX ...
```

Keep `vcfcache annotate` and immutable YAML recipes for custom pipelines,
containers and reproducible advanced workflows. The wrapper should translate
supported invocations into the same internal execution specification rather
than introduce a second annotation engine.

## Annotator adapters and argument handling

- Candidate adapters: VEP, fastVEP and vep-rs. Begin with a narrow, explicitly
  tested subset rather than guessing how arbitrary commands work.
- Each adapter identifies input/output arguments, output format, reference and
  transcript resources, plugins, supplementary databases and semantic options.
- Accept argument ordering supported by the underlying annotator. Preserve
  repeated arguments, plugin order and last-option-wins behaviour; never sort
  arbitrary arguments to construct a fingerprint.
- Treat everything after `--` as annotator arguments, except the adapter's
  deliberate input/output rewriting for the miss-only invocation. Do not
  reconstruct commands through unsafe shell-string interpolation.
- Unsupported options or unresolvable semantics must not silently enable
  caching. Report the limitation and offer an explicit recipe-based route.

## Automatic cache discovery, strict compatibility

Resolve the invocation to an annotation fingerprint. Compatibility must include
the annotator identity/version, assembly, reference and transcript resources,
plugin code and data, supplementary databases and annotation-affecting options.
Only exclude options proven irrelevant to annotation semantics for that adapter.
Do not infer assembly from chromosome naming alone.

Search local caches first. An exact match may be selected automatically; explain
the selected cache and provenance. If several exact matches differ in coverage
or size, make the selection policy visible and allow an explicit override.
Never substitute the nearest-looking recipe or a different annotation version.

Resource identities should be registered and reusable rather than rehashing
large databases on every sample. Pin immutable resources where possible, and
define invalidation/revalidation rules for mutable local files; a path alone is
not a content identity.

## Predictable fallback and resource acquisition

- With no compatible cache, run the original annotator directly and explain why
  caching was unavailable. Preserve its exit status and diagnostics.
- A proposed `--require-cache` option should turn unavailable or incompatible
  caches into a clear error instead of a direct run.
- Large downloads, new cache builds and persistent cache updates require an
  explicit user action or configured policy. Discovery is not permission to
  download or build.
- Base caches and their recipes remain immutable. Any future cohort-learning
  overlay needs its own version, manifest and recorded identity.

## Preserve the output contract

Respect the requested VCF/BCF format and annotation content, with the existing
sorted/indexable VCFcache output guarantee made explicit. Validate complete
records and relevant headers against the corresponding direct invocation.

Do not claim arbitrary commands are safely cacheable. Genotype-, sample- or
cohort-dependent annotations may not be reusable by variant key alone. Auxiliary
reports and statistics produced by annotating only misses are not automatically
equivalent to full-sample reports. Adapters must either reproduce such outputs
correctly or reject/explain the unsupported contract; do not silently emit
partial reports as if they covered the whole sample.

## Other high-value interface improvements

1. **Explain / dry run.** A proposed `--explain` mode should show the detected
   pipeline, resolved resources, fingerprint, selected cache, compatibility
   checks, fallback reason and planned execution without annotating the sample.
2. **One CPU budget.** Expose a top-level thread budget and divide it between
   routing, annotation and compression without oversubscription. More cores
   should be usable where beneficial, not advertised as guaranteed linear
   scaling. Keep temporary-directory and disk requirements discoverable.
3. **Reusable pipeline registration.** Resolve tool and resource identities
   once, save an inspectable pipeline specification, and reuse it across
   samples with explicit invalidation when dependencies change.
4. **Clear progress and provenance.** Keep diagnostics on stderr, preserve
   machine-readable stdout, and record the selected cache, actual command,
   engine, CPU allocation and whether execution used caching or direct fallback.

## Relationship to performance work

The wrapper improves adoption; it does not itself remove startup costs. Share
the fused-routing backend, avoid repeated discovery and unnecessary conversions,
and measure startup on small inputs as well as throughput on WGS. A future
cost-based bypass should use validated measurements and explain its decision,
not silently assume caching always helps.

Before implementation, settle the first adapter and its supported output
contract, compatibility schema, cache-selection policy and fallback behaviour.
Test argument parsing, resource changes, mismatch handling, output equivalence,
exit codes and CPU budgets before extending to additional annotators.
