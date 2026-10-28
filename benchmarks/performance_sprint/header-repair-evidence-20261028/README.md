# Header repair launch evidence

These are launch/validation snapshots, not a claim that all repairs finished.

- Deployed source `787a3081a`, compatible helper 0.7.1.
- Host and Docker: **666 passed, 21 skipped**; deployed GitHub CI passed.
- Gate 2001: **165 complete records, zero mismatches, relevant headers equal,
  no ignored fields**. This unstarted final panel condition was measured once.
- `AUTO_DISPATCH.json`: **35** accepted replacement jobs, 2004–2038; collector
  2039. `FINAL_PARENT_COLLECTOR.json`: 2040 depends on original 2000 and repair
  2039, so the final review covers the complete original 600-slot inventory.
- `CONTROLLER_DISPATCH.json`: metadata-only jobs 2002/2003 were superseded after
  the gate passed. Collection/dispatch ran on the controller rather than waiting
  for an exclusive worker. No successful annotation was repeated.
- `DIRECT_REUSE_PREFLIGHT.json`: three completed original WGS direct controls
  passed full artifact/binding verification; their 17,265.99 seconds (4.8
  node-hours) of completed annotation are preserved. A fourth in-flight control
  is protected by dependency `afterany:1568` on replacement 2007.
- `replacements-a3.json` / `replacements-a4.json`: pending original counterparts
  were cancelled only after replacements were accepted. No active original was
  killed. All original frozen sources, caches and failed attempts remain intact.
- Repair job priority was raised for prompt execution, with the direct-reuse WGS
  retries first. CPU/memory settings and node exclusivity were not changed.
- At 21:23 UTC, WGS retries 2004/2005 were active, all six workers were busy, and
  worker-side 2004 had successfully reused original direct job 1548. Its cached
  sampler was active. Remaining validation is intentionally not predeclared green.

Refresh current progress using `ACTIVE_MANUSCRIPT_CAMPAIGN.md`. The authoritative
checkpoints are the original campaign's higher-numbered a003/a004 attempts.
