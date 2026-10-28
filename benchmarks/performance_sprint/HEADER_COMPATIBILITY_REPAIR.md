# INFO cardinality recovery — 28 October 2026

The user approved a targeted repair after three clinical SnpEff GRCh37 WGS
conditions failed before cached annotation. Their direct controls succeeded and
are retained. The failure was `AC: Integer/. versus Integer/A`: the receiving
input header permits variable-length vectors but the cache declares per-ALT
integers. Rejecting this representable overlay was overly strict.

## Product change

- Preserve the receiving header and all data values. A receiving `Number=.` may
  accept a narrower incoming cardinality **of the same INFO type**.
- Do not narrow a variable incoming vector into a fixed receiving declaration.
  Type conflicts and distinct fixed cardinalities still fail before overlay.
- Python and the optional resident helper implement the same directional rule.
  Native helper **0.7.1** is required when using this broadened compatibility.
  Adaptive routing does not change its per-record algorithm.
- Real bcftools differential tests cover fused/adaptive/native engines with
  zero, partial and complete hits; complete records and relevant headers match
  legacy, including preservation of `AC Number=.`. Negative cases remain errors.

## Campaign recovery

The original 600-condition manifest, old deployments, caches and all failed
attempts are preserved. The potential scope is exactly the 36 clinical SnpEff
GRCh37 conditions. Pending originals were held; active originals were not killed.
Already successful conditions are excluded from replacement.

`repair_headers.py` freezes a separate repair manifest and build, then executes
new higher-numbered attempts in the original checkpoint tree. Each replacement
records its build commit and source-manifest digest. All other conditions retain
their original frozen build. The delta changes only header compatibility checks,
not annotation commands, input normalization, routing algorithms or CPU budgets.

Completed direct controls are reused only after verifying the exact source task
(input, assay, parameters, CPU and resource binding), original command digest,
successful direct resource metrics, affinity, checkpoint and output/index hashes.
An in-flight original control is allowed to finish; its successor depends on
that job and freezes the completed artifacts before cached timing. A changed
condition or unsuccessful control cannot be reused silently.

One actual final panel condition gates automatic release of the remaining
repairs. New accepted job IDs are registered in the original submission history;
only corresponding pending predecessors are cancelled. A final collector waits
for both the original campaign collector and the repair collector. No successful
sample/tool condition is repeated, and no direct annotation is repeated solely
because the cached header guard failed.

Launch with a fresh receipt and the compatible helper:

```sh
python -m benchmarks.performance_sprint.launch --helper COMPATIBLE_071_HELPER \
  --receipt UNIQUE_REPAIR_RECEIPT.json \
  --repair-parent /results/manuscript-final/20261027T191111Z-fe197f34d/campaign
```

Submission receipts and the original campaign's live checkpoints, not this
document, establish whether recovery is running or complete.
