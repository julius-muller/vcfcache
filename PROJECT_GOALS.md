# Project goals

VCFcache makes repeated variant annotation faster by reusing annotations that
were produced by exactly the same annotation recipe. Its primary contract is
that caching changes runtime, not results.

## Product goals

1. **Preserve annotation results.** Given the same normalized variants,
   annotation recipe, reference data, and tool version, cached and direct runs
   must produce annotation-equivalent output. Performance claims are invalid
   until this has been checked across every output record.
2. **Remain useful as annotators become faster.** VCFcache must reduce total
   end-to-end runtime substantially at realistic WGS hit rates, including for
   native tools such as fastVEP and future Ensembl annotation engines. Faster
   annotation must not make cache lookup, preprocessing, or merging the
   dominant cost.
3. **Keep the public workflow tool-agnostic.** Command recipes remain the
   general integration boundary. Optimized native adapters are appropriate
   when a command-based boundary cannot provide meaningful acceleration, but
   they must not make existing tools or published caches second-class users.
4. **Make performance claims comparable and intuitive.** A cached benchmark
   and its direct baseline must use the same input, preprocessing policy,
   annotator, databases, output contract, hardware, and resource limits.
5. **Support both public and cohort caches.** Users should be able to understand
   the benefit for a new sample from its hit rate, whether the cache originates
   from a published resource or a previous cohort.
6. **Spend compute on independent evidence.** Publication campaigns execute
   each unique sample, annotator, and cache condition exactly once. Adapter
   smoke tests use a small separate input; full-genome warm-ups and technical
   reruns are not part of the design.

## Fast-annotator performance target

The initial target for a matched whole-genome pipeline is at least 2x speedup
at an 80% cache-hit rate and at least 3x at 90%, without changing output. A 50%
hit-rate run should not be more than 10% slower than direct annotation. These
are product targets, not guarantees for every annotator or storage system.

The first core fastVEP experiment and a follow-up dense supplementary-database
stress pilot are exploratory engineering results. The dense configuration met
the target at 90% hits (3.60x) with complete output equality, while the core
configuration exposed a lower ceiling. Both identify the same fixed-cost
passes: lookup output, miss selection, and final output construction. A
streaming engine may subsequently be developed inside VCFcache or as a
separate wrapper, upstream contribution, or fork around a native annotator.

The publication follow-up treats fastVEP as a peer annotator, not a special
pilot: it uses the exact 52 held-out KPGP, SGDP, and PGP genomes, the same
bundled gnomAD AF thresholds and disjoint three-genome cohort blueprints, the
same balanced four-condition task layout, and the same output-equality gate as
VEP. The blueprint membership is shared; the annotations stored in each cache
are rebuilt with an assembly-matched fastVEP recipe.

## Ensembl platform transition

Ensembl 116 is the final release on the current Ensembl platforms. This is not
treated here as evidence that VEP 116 is the final VEP release. VCFcache will
track the next Ensembl annotation platform and evaluate it using the same
correctness and matched-pipeline performance gates when it becomes available.
