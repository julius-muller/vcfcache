# Author-controlled submission items

These items are intentionally not inferred from the repository.

## Required before bioRxiv submission

- [ ] Add a representative non-perfect-hit WGS stage profile with measured
      preprocessing, lookup, miss selection, annotation, merge and write/index
      durations; the existing perfect-hit diagnostics are not an additive profile.
- [ ] Audit assembly, contig and reference-normalization provenance for each
      input/cache pairing; do not infer left alignment from `bcftools norm -m-`.

- [ ] Final title.
- [ ] Complete author order, affiliations, email addresses, and ORCIDs.
- [ ] Corresponding author designation.
- [ ] Institutional wording for analysis of openly accessible,
      de-identified public human genomic data.
- [ ] Funding bodies, grant numbers, and funder roles.
- [ ] CRediT author contributions.
- [ ] Competing-interest declaration.
- [ ] Acknowledgements and permission from named contributors.
- [ ] Final VCFcache release version, Git commit, container digest, and archive
      DOI used to support the manuscript.
- [ ] Permanent DOI for benchmark source tables, scripts, manifests, and
      sufficient test input/output to reproduce the examples.
- [ ] BioRxiv license selection.
- [ ] Name of the author who directed the disclosed AI-assisted drafting.

## Required or strongly recommended for GigaScience

- [ ] Register or confirm the VCFcache bio.tools identifier.
- [ ] Register or confirm an RRID, if applicable.
- [ ] Freeze an archival code snapshot in Zenodo, Software Heritage, or
      GigaDB.
- [ ] Decide whether to deposit a compact executable example in Code Ocean.
- [ ] Confirm that test inputs and expected outputs can be freely downloaded
      by editors and reviewers.
- [ ] Finalize the GigaScience cover letter.
- [ ] Recheck GigaScience file naming and submission requirements immediately
      before submission.
