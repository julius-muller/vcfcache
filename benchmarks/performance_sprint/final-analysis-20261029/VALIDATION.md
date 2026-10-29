# Analysis validation

29 October 2026. Verification concerns this consolidation, not a new annotator
or release validation campaign.

- Full local Python suite: **680 passed, 21 skipped**. Fourteen new analysis
  tests cover condition/alias counts, retained attempts, exclusions, phase
  identities, archive hash, training-seed leakage, strict repayment thresholds,
  step counters, immutable snapshot protection and fail-closed comparison audits.
- Standalone R checks passed: singleton-stratum bootstrapping, no alias/sample
  inflation, measured phase fractions summing to one, and refusal to summarize
  repeated donors as independent observations.
- New Python sources pass Ruff and Black formatting checks.
- R rendering completed without warnings/errors. Plotting uses ggplot2 and
  patchwork, with base Cairo PDF/PNG and SVG devices. Dependencies/session are
  recorded in `source/R_sessionInfo.txt`.
- Representative full panels were visually inspected: assay performance,
  WGS strategies, phase/headroom, synthetic load, CPU/memory and recorded build
  amortization. All composite panels also have individual exports.
- Preferred `figures-v2/` contains 11 complete graphics and eight individually
  exported lettered panels, each in PDF/SVG/600-dpi PNG. Original `figures/`
  exports and all pre-existing repository figures remain intact.
- Source archive SHA256 matches `source/SNAPSHOT.json`; individual archived
  metadata hashes were rechecked. Cohort seed IDs are disjoint from their
  strategy test donors. gnomAD donor-disjointness was not established.
- Build amortization uses provenance matching each executed cache payload hash;
  stale strategy-manifest build durations remain alongside the selected costs.
- No failed or timing-only output was promoted by this analysis. Counterfactual
  headroom and amortization curves are explicitly distinguished from observed
  annotation measurements.

The archive intentionally excludes large BCFs and dense sampler streams.
Their original retained paths and frozen hashes remain available; this local
validation did not rerun a complete semantic scan of every remote BCF.
