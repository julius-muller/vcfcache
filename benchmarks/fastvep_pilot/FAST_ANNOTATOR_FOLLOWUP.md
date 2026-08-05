# Possible native fast-annotator follow-up

This is a decision record for work after the exploratory fastVEP pilot. It does
not commit the current VCFcache release to a Rust rewrite.

## Decision trigger

A follow-up is justified when cached fastVEP output is correct but the current
multi-pass VCFcache workflow misses either product target: 2x speedup at 80%
hits or 3x at 90% hits. The pilot's 100%-hit stage profile determines where the
fixed cost originates.

## Alternatives

### Native backend within VCFcache

Keep the existing CLI, recipes, and cache archives while adding a native engine
that performs lookup and output assembly with fewer full-file passes. Prefer
this when BCF lookup, intermediate writes, and merging dominate and the same
engine can accelerate several annotators.

### Standalone wrapper around fastVEP

Build a separate Rust application that reads each variant once, performs the
VCFcache lookup, sends only misses through fastVEP, and writes one ordered
output stream. Prefer this when a single-process pipeline is required but a
clean boundary around upstream fastVEP remains possible.

### fastVEP contribution or fork

Integrate annotation-cache lookup directly with fastVEP's internal Rust crates.
Prefer an upstream contribution; use a fork only if the required interfaces are
not accepted or exposed. fastVEP declares Apache-2.0 licensing and currently
separates cache, I/O, consequence, HGVS, and annotation crates, so the approach
is technically plausible. It also couples VCFcache to a young and changing API.

## Selection rule

- Choose the generic VCFcache backend when generic BCF passes dominate and a
  prototype can meet the speed target.
- Choose the wrapper when serialization and process boundaries dominate but
  fastVEP can remain an upstream dependency.
- Choose contribution/fork work only when direct access to fastVEP's in-memory
  transcript and output structures is necessary to meet the target.
- If fastVEP has material core-consequence gaps, retain the native-engine design
  for future annotators but do not present fastVEP as a current VEP replacement.

Any follow-up must preserve existing Zenodo cache archives or derive optional
local sidecar indexes from them. It must prove output equality before changing
the default engine.

