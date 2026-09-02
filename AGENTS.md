# Repository Guidelines (untracked helper)

## Structure
- Code: `vcfcache/` (CLI, database orchestrators, utils, manifest helpers).
- Assets: `resources/`, Docker files in `docker/`, docs in `README.md`, `WIKI.md`, manifest in `public_caches.yaml`.
- Tests/fixtures: `tests/`, configs in `tests/config/`.

## Build/Test
- Dev env: `uv venv .venv && source .venv/bin/activate && uv pip install -e ".[dev]"`.
- Tests: `python -m pytest` (network/Zenodo integration gated by env).
- Lint/format/type: `ruff check --fix vcfcache`, `black vcfcache`, `mypy vcfcache`.
- Docker: multi-stage `docker/Dockerfile.vcfcache`; `--target test` installs dev deps and can run pytest; `--target final` is the lean runtime (`ghcr.io/julius-muller/vcfcache:<version>`).

## Manifest & Aliases
- Aliases: blueprint `bp-<genome>-<source>-<release>-<filt>`; cache `cache-<genome>-<source>-<release>-<filt>-<tool>-<tool_version>-<preset>`.
- Manifest entries include `version` (vcfcache used), no image tags. Fields: alias, type (blueprint|cache), genome/source/release/filt, doi, updated_at, md5; caches add tool/tool_version/preset and annotation_yaml_md5.
- CLI resolves aliases via `public_caches.yaml`; caches/blueprints live on Zenodo.

## Workflow tips
- Keep `annotation.yaml` immutable per cache; store `params.yaml` snapshots with caches.
- Record new aliases in `public_caches.yaml` when publishing.

## Docs
- `WIKI.md` is canonical for alias/manifest schema, Zenodo pull/push, Docker usage, cache layout. README mirrors high-level usage.

## Git workflow
- After completing and verifying requested repository changes, commit them and push the current branch unless the user explicitly asks not to.
- Author commits as `Julius Müller <julius.mueller@dkfz-heidelberg.de>` rather than with the VM's default `ubuntu` identity.
- Inspect the worktree before staging, preserve unrelated user changes, and describe the committed scope in the final response.
