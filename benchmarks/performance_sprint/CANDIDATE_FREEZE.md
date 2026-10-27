# Candidate selected after decision A — 27 October 2026

The user selected **A: stop optimizing and execute the full manuscript
benchmark**. No additional generic approval is needed. The 36-cell optimization
screen is closed; the following are publication measurements, not an expanded
engineering screen.

## Frozen implementation choices

- Primary caching engine: **adaptive**, across tools and assays. Native remains
  opt-in; its existing WES evidence is not silently pooled into the primary arm.
- Enable the bounded early-fallback guard for native requests that would
  otherwise fill an over-budget resident buffer. Its measured benefit is limited
  to avoiding failed resident attempts; it is not a universal adaptive speedup.
- Leave indexed-miss, smaller-window, trusted-hit and bounded-router experiments
  disabled. Do not promote sparse-input wins into a universal routing default.
- Keep the pinned annotator binaries/recipes and original matched CPU settings.
  No upstream vep-rs upgrade, bug investigation or new cache format.

## Coverage and baseline reuse

Materialize all **540 historical cached rows**, deduplicating only byte-/command-
identical execution/reference bindings and retaining every source ID as an alias.
Expected overlap is twelve VEP WGS rows also used in the assay-extension plot.
Add **72 new paired vep-rs conditions**: twelve GRCh38 genomes, three assays,
core and clinical recipes. Each new pair measures its own direct baseline once;
historical VEP/fastVEP/SnpEff direct measurements are not routinely repeated.

Approvals bind SHA-audited input/cache/recipe/resource files, exact runtime CPU
settings and retained reference evidence. A successful SHA audit can be reused
only with unchanged file size and nanosecond mtime; a metadata change triggers a
full SHA check. Current file fingerprints do not imply contemporaneous historic
OS/page-cache state. Measurements retain their baseline provenance and timing
scope rather than pretending all direct measurements are new.

- Retained outputs: strict new-versus-cached comparison; direct comparator stays
  annotator-specific. VEP's documented CSQ exception is permitted only alongside
  complete non-CSQ equality and equal relevant headers; raw differences persist.
- C2 digest-only rows: reproduce their original exact annotation/key digest and
  label that limited validation scope. Preserve the known SnpEff exclusions.
- Synthetic loads: delay emits no annotation. Bind historical timing/semantic
  metadata and compare against the same sample's retained VEP annotation recipe.
  Do not claim that missing per-delay direct BCFs were retained or revalidated.
- New vep-rs core: strict fresh-pair comparison. Clinical CSQ discrepancies remain
  timing-only, never validated speedup claims; differences outside CSQ fail.

## Autonomous execution after submission

1. Immutable deployment plus an eight-CPU preparation job freezes all resources,
   builds the complete task document and checks baseline bindings. No annotation
   begins if preparation fails.
2. Submit six representative **actual primary conditions**, not extra repeats,
   across all six eligible workers. Each timed condition has its own Slurm step.
3. An afterany collector persists their states. A dependent dispatcher releases
   the remaining matrix only if those conditions pass their declared checks.
   Otherwise it writes `DISPATCH_BLOCKED.json`; this is a concrete blocker.
4. Remaining jobs queue automatically, longest historical cached conditions
   first, and refill freed nodes. A final afterany collector records the full
   status, including failures. Successful conditions are never silently rerun.

The launch receipt records preparation acceptance; `START.json` records the
first six jobs and continuation; `AUTO_DISPATCH.json` records bulk submission.
These files, not this document, prove which stages have actually executed.

No user decision is expected during a healthy measurement campaign. The next
decision is after the collected results/updated figures, or earlier if a concrete
validation/resource failure requires a change outside the approved plan.
