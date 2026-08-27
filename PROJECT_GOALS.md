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

The initial target for a whole-genome annotation pipeline is at least 2x speedup
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

The publication follow-up evaluates VEP, fastVEP and SnpEff with realistic
native configurations. The primary estimand is always the speed-up from direct
to cached annotation within the same tool, recipe, input and hardware. Absolute
runtime is not used to rank the annotators: their transcript resources,
consequence models, HGVS implementations, supplementary fields and output
schemas are not identical. VEP `--everything`, fastVEP core and supplementary
database configurations, SnpEff core and clinical configurations, and the
synthetic-delay sweep deliberately span the pipeline-cost range users encounter.
The same bundled gnomAD AF>=1% blueprint membership is used where applicable,
but annotations stored in each cache are always produced by that cache's own
annotator recipe.

The inferential units are independent genomes, not repeated technical runs.
Each unique sample/tool/assay/condition is timed once; external semantic
comparison is performed afterward and outside the timed cell.

## Ensembl platform transition

Ensembl 116 is the final release on the current Ensembl platforms. This is not
treated here as evidence that VEP 116 is the final VEP release. VCFcache will
track the next Ensembl annotation platform and evaluate it using the same
correctness and matched-pipeline performance gates when it becomes available.
