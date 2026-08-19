# Independent matched assay extension

This campaign supports v3 Figure 2 without reusing the historical
1000 Genomes-derived Panel/WES inputs. It selects six KPGP and six SGDP GRCh38
genomes from the independent external evaluation set and derives two matched
inputs from every genome:

- ACMG SF v3.3 MANE Select coding regions with 20 bp padding (`Panel`);
- Twist Human Core Exome capture regions with 125 bp padding (`WES`).

Each input is measured once with VEP 115.2 and fastVEP 0.3.0, directly and
through the gnomAD AF ≥1% cache strategy. There are no technical repeats.
The same 12 genomes' already-frozen WGS observations are joined during
collection, producing a complete 2 annotators x 3 assays x 12 genomes table.

The public VEP cache is the production Zenodo artifact. Because Zenodo does
not distribute a fastVEP-specific annotation cache, the fastVEP cache is the
locally annotated counterpart built from the frozen bundled blueprint; its
provenance remains distinct from that of a downloadable cache.

## Reproduction components

- `../prepare_external_assay_v3.py`: deterministic selection and manifests;
- `../slurm_prepare_external_assay_v3.sh`: atomic interval subsetting;
- `../run_external_assay_task.py`: paired conditions and annotator-specific
  semantic comparison;
- `../slurm_external_assay_v3.sh`: separate Slurm job steps per condition;
- `../collect_external_assay_v3.py`: strict 72-row matrix collection.

The frozen campaign ID is `external-assay-v3-independent-af001-warm`. Generated
manifests are retained below the same-named subdirectory. After completion,
machine-readable publication data are written under `source_data/`.

Because Panel/WES timed cells are short relative to annotator and filesystem
cold-start effects, every task first processes the same separate KPGP-00319
assay subset through the cached path. This warm-up is excluded from reported
wall time and is not an evaluation-sample replicate. It places the annotator,
transcript assets and cache lookup path in a consistent warm state before the
single measured direct and cached conditions; their order remains balanced.

## Interpretation

The assay extension is deliberately a use-case boundary experiment. A tiny
panel can be faster to annotate directly because normalization, cache lookup,
merge and indexing are fixed per-sample costs. That outcome is retained rather
than hidden. Figure 3 separately establishes how the balance changes for more
expensive annotation recipes.
