# Bounded Rust WES engineering pilot

This is an isolated experimental branch, not a product-default change or a new
publication campaign. The two targets are **fastVEP clinical WES** and **vep-rs
core WES**, using the existing immutable GRCh38 recipes and caches. A real cached
speedup is required; silently choosing direct annotation does not count.

## Scope and acceptance

- Twelve independent WES inputs per recipe, copied from the frozen completed
  manuscript matrix. No new assays, annotators, annotation recipes, or formats.
- Every timed CLI invocation gets a separate eight-CPU Slurm step on an
  exclusive worker. Preparation and strict complete-record/header comparison
  are outside timing. Statistics use the product's light mode.
- Primary target: cold **application-start** break-even on every input;
  speedup of 2× is a stretch target. Filesystem caches are **not** flushed.
- No ignored INFO/FORMAT/CSQ fields, no truncated outputs, and no native fallback
  may pass. The direct output must also match its frozen retained reference.
- Session experiments must disclose loading, aggregate memory, and amortization
  separately; their warm latency must never be substituted for cold CLI latency.
- Three-hour wall-clock cap from 2026-11-03 09:44:15 UTC. Whole-cache loading was
  abandoned once capacity probes showed it could not fit the available budget.

## Configurations

| Name | What it measures |
| --- | --- |
| Direct | Same recipe through `annotate --uncached`; complete CLI baseline |
| Rust-adaptive | Existing sparse fused-routing implementation |
| Resident input | Existing native implementation retains the **sample**, not the whole cache; default four lookup workers at eight CPUs |
| Coalesced input | Experimental resident-input lookup merges windows whose actual CSI virtual-offset spans overlap, with an 8-MiB compressed-span guard |
| Eight lookup workers | Same resident-input algorithm with `--engine-workers 8`; zero extra cache-decoder pool threads under the same eight-CPU budget |
| Whole-cache session | Experimental lazy standard-BCF records retained across inputs; **capacity test failed**, so no validated warm-throughput claim |

The coalescing switch is experimental:
`VCFCACHE_EXPERIMENT_COALESCE_CHUNKS=1` (helper >=0.8.0). Normal native defaults
are unchanged. The pilot-only helper commands `resident-batch` and `cache-probe`
use a 20-GiB sampled **helper-RSS guard**, not a 20-GiB aggregate hard limit;
Slurm jobs have a separate 28,000-MB memory allocation.

The controller executes the recorded Bash script with `set -euo pipefail` so a
failed upstream pipeline cannot be reported as successful. This safety fix is
independent of the experimental lookup choices.

## Freeze, submit, and inspect later

Run `prepare` on the Slurm controller with the frozen publication dataset
available. Deployment is an immutable checkout; the helper path must be frozen
as well. Each new output root must not exist already.

```sh
python -m benchmarks.wes_pilot.prepare /mnt/data/slurm-results/wes-pilot/NEW-ROOT \
  --deployment /results/deployment/FROZEN-CHECKOUT --commit FULL-COMMIT \
  --binary /results/deployment/FROZEN-BINARY/vcfcache-engine
```

`--no-submit` creates the immutable manifest and READY hash without starting
jobs. `--selection remaining --balance-order` schedules the other nine inputs
per recipe after the deterministic three-input screen. `--exclude-ids` avoids
rerunning already completed confirmations. For an additional configuration,
`--reference-screen /results/wes-pilot/COMPLETED-ROOT` reuses its measured direct
baseline, while output comparison uses the SHA-bound retained direct output
that the fresh direct already matched. It is not a simultaneous direct repeat.

Manifests, submitted commands/job IDs, per-condition CPU/memory measurements,
states, strict comparisons, output hashes, and errors are durable. Rejected
attempts remain visible; the collector never picks the fastest technical repeat.

## Collect and plot

The dated results directory contains a compact terminal-evidence archive and
SHA-256 inventory. Extract its `source/` member in that directory, then:

```sh
python -m benchmarks.wes_pilot.collect benchmarks/wes_pilot/results-20261103
Rscript --vanilla benchmarks/wes_pilot/render.R benchmarks/wes_pilot/results-20261103 NEW-FIGURES
```

Plotting is exclusively R/ggplot2. The renderer refuses an existing figure
directory and exports vector PDF/SVG plus 600-dpi PNG, including individual
panel files. Numeric observations and R tibbles are available as TSV/RDS.

See [RESULTS.md](results-20261103/RESULTS.md) for the decision, caveats, failures,
and exact results. These are engineering measurements, not a replacement for
the manuscript campaign or evidence of linear CPU scaling.
