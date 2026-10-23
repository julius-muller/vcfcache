# Frozen pilot source snapshots

These copies identify the sources actually deployed, rather than substituting
the final feature-branch tree for earlier pilot code. Exclude them from normal
development/linting; edit the maintained sources in the repository root.

- A: cache-derived size/span controls, jobs 1257–1259.
- B: real WGS and experimental store, jobs 1261–1262.
- C: untimed cache-payload read trace, job 1265.
- F: final real panel/WES validation, job 1268.

Deployment directories are `/results/deployment/vcfcache-adaptive-pilot-20261023X`
on workers, with lowercase X. Executables/build trees and environments are
excluded; run manifests record executable/source hashes. Early native builds
retained the development 0.1.0 label; the final helper requires version 0.2.0.
Original fused baseline is commit c878d86c3, with native executable SHA256
`c48c5f8633b82c3902296fb1e777649b22ac0e5b3edf722871213015d38979ad`.
