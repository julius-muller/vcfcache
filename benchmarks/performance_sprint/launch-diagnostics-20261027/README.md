# Manuscript launch validation, 27 October 2026

These are technical launch diagnostics, not additional manuscript replicates.
Large annotated BCFs remain on the cluster under
`/results/manuscript-final/20261027T135803Z-b4b43512c/`.

- Overlay audit 1386 reused the failed checks' annotations, without running an
  annotator. The adaptive output and reconstructed legacy INFO overlay compare
  equal for both clinical examples. This isolates the issue from routing.
- Independent direct controls 1387 used the original bare C2 wrappers. fastVEP
  KPGP-00128 panel passed complete-record/header comparison. Its new direct
  projection digest matched the new cached digest, not the archived projection.
  No historical digest was silently replaced or treated as a successful audit.
- SnpEff KPGP-00139 panel initially differed at 239/242 records, largely because
  SnpSift adds reference identifiers to the VCF ID column by default. An untimed
  diagnostic removing ID on both sides still found **one** INFO difference:
  the accented clinical disease name at chr2:47476361 (CLNDN). This diagnostic
  did not authorize ignoring either field in the final comparison.
- The existing cache contains `Muir-Torr??_syndrome`, whereas the independent
  UTF-8 direct control produced `Muir-TorrÃ©_syndrome`. Neither establishes Unicode
  fidelity to the source database. Java 21 defaults to UTF-8 even with LC_ALL=C.
  The compatibility recipe therefore explicitly pins ASCII text encoding and
  SnpSift `-noId` in **both** arms. This preserves the existing cache's text
  representation; it is not a Unicode-correctness claim or an annotator repair.
- Corrected direct control 1394 passed complete-record/header comparison on all
  242 records, with **zero ignored fields**, against the retained cached output.
  This control changed configuration and is not a technical repeat for averaging
  or choosing the fastest runtime. Its 122.96 s is diagnostic evidence, not a
  replacement primary campaign timing. The final matched recipe is measured once.

Worker-local fastVEP caches were independently generated. A two-node probe found
7,967,504 records on both nodes but different compressed BCF sizes and hashes;
the provenance JSON differs on all six workers. The recovery freezes canonical
payload, index, recipe, parameters and provenance together. Drifted copies are
installed into private node-local directories outside timing. Existing caches
are never overwritten; compatible original replicas are reused.

The incomplete preparation 1395 was explicitly cancelled before annotation after
the provenance-replica issue was confirmed. Its receipt and deployment remain
archived. Preparation 1397 uses the corrected canonical-provenance staging.
