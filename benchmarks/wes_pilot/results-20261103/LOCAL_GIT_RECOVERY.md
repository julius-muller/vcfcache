# Local packaging incident

During final packaging on 3 November 2026, recently readable Git objects became
unavailable in the shared local object store. The cause is not established.
Benchmark working files, Slurm evidence, and immutable controller deployments
were unaffected. A backup of the damaged Git metadata was made under
`/tmp/vcfcache-git-recovery.KMtrPV/damaged-git` before recovery; no reset, checkout
replacement, pruning, or history rewrite was performed.

Published objects were restored additively from a fresh mirror of the existing
origin. Missing unpublished pilot blobs/trees were recovered from the immutable
deployed source and the recorded original patches. Commits were reconstructed
with their original trees, parents, authors, messages, and timestamps. Every
recovered commit matches its **original exact SHA**; this is not a substituted
history or a retroactively changed benchmark source version:

```
f61672e3ef0a19feeaf09e7cf31b4506f5c40283
0a0ef66d86282fe340e79294e05466d24a38e69a
5565de6a516d75c71c9d30a5619c8991a4c08c0f
34b87036fe8f51cabf6baf53b57e88cb8c3614e8
55dfe0395102b44c5ec6a2f22f613123a3036d35
d14a5012487f154ca12a6f39831fe256d14a243c
ac67f46afa10976aaeb885bd8c086363d1386fa0
0f0a9eee944c0a1724c6ba62248c8b25a3d6e05a
59b49933fad4bf334be49cce228acd5187fdce6f
```

The original `feature/python-adaptive` worktree remained clean and unswitched.
Recovery copies are intentionally retained, not deleted. No unrelated
synchronization service or user process was stopped or reconfigured.

All current references and their history pass
`git fsck --full --strict --no-reflogs --no-dangling`. Including historical
reflogs additionally reports an unavailable old Python-adaptive commit,
`5228090d131bce30c721906875f6bdd155dd9a23`. It is not reachable from any current
reference, and origin cannot supply it by SHA. Those unrelated reflog entries
were deliberately left intact; recovery does not claim to have restored that
old, unreachable commit or established why objects disappeared.
