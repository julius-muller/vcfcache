# Core, annotation-rich pipelines and model sensitivity

The final outputs are in **figures-v3/**; earlier local drafts are not publication
assets. The original final-20261103 dataset and its figures are unchanged.

## Is core WGS fundamentally hopeless?

No. On the same 12 samples, the frozen **Rust-adaptive** campaign gives median
vep-rs core WGS speedup **1.32×**, compared with **0.89×** for Python-adaptive.
All 12 Rust core WGS conditions are faster than their own direct baseline
(range 1.09–1.49×). Core WES loses in both implementations (0.46× median).
This is historical evidence of headroom, not a fresh randomized engine comparison.
The direct denominators, sample and cache hashes are identical across the two sets;
execution dates, nodes and OS/filesystem-cache state can still affect the difference.

For Python core WGS, the median of per-sample runtime reductions needed merely
to reach parity is **11.2%**. That is a modest target, not a demonstrated improvement.
The current logs also show avoidable-looking data movement rather than a lookup-only
problem. Across 12 core WGS conditions the median logged command-wall contributions
are approximately 21.2 s normalization, 16.1 s cache overlay, 18.2 s miss annotation,
8.6 s merge overlay, 12.5 s gather/index and 17.2 s final write/index. These medians
must not be added as if they describe one sample. Shard wall times use the maximum,
not the sum of concurrently executed commands. Residual/coordination time is retained
in `tables/core_WGS_log_headroom.tsv`. This is a **log-derived critical-path
approximation**, not an instrumented phase profile.

### Recommended bounded engineering follow-up

1. Try filtering and final compression within each merge shard, then gather final
   shards into the requested output: avoid the extra full-size merged intermediate
   followed by another complete read/filter/write. Preserve sorting, header union,
   allele pairing, output indexing and strict equivalence. Output/compression cost
   cannot disappear; it must be paid once at an earlier stage. A zero-cost deletion
   of the final logged stage is only an optimistic bound (median 1.06×), not a forecast.
2. Test two/four shards at the same eight-CPU allocation rather than assuming that
   more concurrent bcftools processes are always better. bcftools compression threads
   are additional to the main process; account for every concurrently live pipeline
   process and keep CPU controls matched.
3. A trusted-normalized-input fast path can improve latency, but offer the same
   normalization contract to **both direct and cached arms**. The benchmark runner
   already normalizes its prepared input, while both current CLI paths normalize
   again. Removing that cost only from cached execution would not be a fair revised
   comparison. Such a baseline change requires a matched direct calibration.

Prioritize the first two small pilots before a deeper rewrite. Retaining Rust as
an opt-in is already a measured alternative; this work does not silently change
the default or schedule another full campaign.

## Where are the important heavy vep-rs results?

They are now prominent in **provisional/Figure4_veprs_core_and_annotation_rich**,
with individual elapsed times and ratios. The clinical recipe is a real annotation-
rich recipe (gnomAD + ClinVar), not a synthetic sleep experiment:

| Recipe / assay | Python-adaptive | Rust-adaptive | Equivalence |
|---|---:|---:|---|
| Core WES | 0.46× | 0.46× | Accepted |
| Core WGS | 0.89× | 1.32× | Accepted |
| Clinical WES | 1.66× | 1.93× | **Timing only** |
| Clinical WGS | 2.72× | 3.18× | **Timing only** |

The clinical timings were not missing from the data: the 24 WES/WGS conditions were
excluded from the validated headline because CSQ differs, including examples where
the cached record lacks ClinVar-derived values present in direct annotation.
Non-CSQ records and relevant headers pass, but a small number of disagreements
still matters clinically. Do not ignore those fields or call the timings validated
speedups. Next resolve the cache/recipe/annotator subset-consistency cause, or use
an annotation-rich recipe whose complete-record equivalence genuinely passes.
The 12 clinical panel conditions pass and are shown as accepted points.

## Figure captions

**Figure 4: vep-rs core and annotation-rich recipe timings.** A, Cached versus direct
elapsed time for each of 12 GRCh38 samples per recipe/assay, with KPGP/SGDP colours.
The diagonal marks equal elapsed time. B, Per-sample direct/cached ratios, with a
dashed parity line. Filled points passed the recorded equivalence comparison;
hollow clinical WES/WGS points have CSQ differences and are **not validated
speedup evidence**. Direct controls are retained measurements at the same eight-CPU
allocation. Each point is a sample, not a technical replicate. CSQ, consequence
annotation field; WES, whole-exome sequencing; WGS, whole-genome sequencing.

**Figure 5: Hypothetical sensitivity of caching benefit to runtime parameters.**
A, Serial annotation cost per variant at 50%, 80% and 93% cache hits for WES-sized
(82,714 records) and WGS-sized (5 million records) inputs. B, Startup cost paid once
per annotator invocation in **both** direct and cached arms whenever misses remain.
C, Configured CPU count under equal ideal scaling, a poorly scaling cache path or
a poorly scaling annotator. Poor scaling means an Amdahl serial fraction of 0.2;
otherwise the serial fraction is zero. D, Cache-path versus annotation cost per
variant; the black contour marks break-even. All panels are **illustrative model
calculations**, not fitted measurements, predictions calibrated to a named tool,
or a new CPU-scaling experiment. Default assumptions (unless varied): 80% hits,
five million input variants, 100 μs serial annotation cost per variant, 20 μs serial
cache-path cost per input variant, 5 s annotator startup, 1 s additional cache setup,
eight CPUs and ideal scaling. Full-input and miss annotation share the same assumed
mean per-variant cost and scaling; real workloads need not satisfy that assumption.
Cache building/download and multi-sample batching are not modeled.

For effective processing scale \(g_i(p)=1/[\sigma_i+(1-\sigma_i)/p]\), where
\(\eta_i(p)=g_i(p)/p\), the model is

\[
T_{\mathrm{direct}}=S_a+\frac{Na}{g_a(p)},\qquad
T_{\mathrm{cached}}=S_c+\frac{Nc}{g_c(p)}+
\mathbf{1}_{f<1}\left[S_a+\frac{(1-f)Na}{g_a(p)}\right].
\]

Here \(N\) is input size; \(f\) is the cache-hit fraction; \(a,c\) are serial
per-variant annotation and cache-path costs; \(S_a,S_c\) are annotator and additional
cache startup; \(p\) is configured CPUs. If misses remain and both invocations pay
the same startup, \(S_a\) cancels from **absolute time saved**, but increasing it
drives the **ratio** towards one. Startup does not independently change break-even
in this simplified matched-startup model. Unequal annotation/cache scaling can
either narrow or widen the caching margin as CPUs increase.

All plots are R/ggplot2, with individual panels, vector PDF/SVG and 600-dpi PNG.
`tables/interpretation_data.rds` contains measured and model tibbles separately.
`source/core_WGS_logs.tar.gz` and per-log SHA256 values retain the command-wall
evidence. See the matching Rust-adaptive figure set in
`../../rust_adaptive_results/final-20261103/` and its provisional Figure 4 in
`../../rust_adaptive_results/interpretation-20261103/`.
