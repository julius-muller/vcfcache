# CLAUDE.md (untracked helper)

Working notes for the current vcfcache layout.

## Architecture
- Pure Python
- Requires `bcftools >= 1.20` (Docker image installs bcftools; for local installs, provide your own).
- Caches/blueprints are portable tarballs; public aliases resolve via `vcfcache/public_caches.yaml`.
- Aliases: blueprint `bp-<genome>-<source>-<release>-<filt>`; cache `cache-<genome>-<source>-<release>-<filt>-<tool>-<tool_version>-<preset>`.
- Genome-agnostic: works for human and any model organism as long as reference/contig naming is compatible across inputs and your annotation pipeline.

## Dev setup
```bash
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
python -m pytest tests -q   # optional
```

## Docker
- Multi-stage `docker/Dockerfile.vcfcache`:
  - `--target test` installs dev deps; use to run pytest in CI.
  - `--target final` is the lean runtime (`ghcr.io/julius-muller/vcfcache:<version>`).
- Runtime image includes `bcftools` and `tabix`.

## Manifest & Zenodo
- Manifest fields: alias, type (blueprint|cache), version, genome/source/release/filt, doi, updated_at, md5; caches also carry tool/tool_version/preset and annotation_yaml_md5.
- Pull via alias: `vcfcache annotate -a <alias>` (auto-download to `~/.cache/vcfcache/caches/<alias>`).
- Pull by DOI: `vcfcache blueprint-init --doi <doi> -o <dir>` (blueprint) or `vcfcache cache-build --doi <doi>` (cache).
- Push/share: `vcfcache push --cache-dir <dir>` (needs `ZENODO_TOKEN`; use `--test` + `ZENODO_SANDBOX_TOKEN` for sandbox).
- List/inspect local items: `vcfcache list caches --local`, `vcfcache list blueprints --local`, `vcfcache list --inspect <path>`.

## Build/test commands
- Lint/format/type: `ruff check --fix vcfcache`, `black vcfcache`, `mypy vcfcache`.
- Run tests (dev): `python -m pytest tests -q`.
- Docker test stage: `docker build --target test -f docker/Dockerfile.vcfcache -t vcfcache:test .` then `docker run --rm vcfcache:test /opt/venv/bin/python -m pytest tests -q`.

## bcftools Configuration
- Minimum version: 1.20 (for `--write-index` flag)
- Version checked at startup: `check_bcftools_installed()` in `utils/validation.py`
- Override system bcftools: Set `VCFCACHE_BCFTOOLS` environment variable
  ```bash
  export VCFCACHE_BCFTOOLS=/path/to/bcftools-1.22
  vcfcache annotate ...
  ```
- Search order: `VCFCACHE_BCFTOOLS` env var → system PATH → error
- For annotation commands: Use `bcftools_cmd` in params.yaml (separate from internal checks)

## Notes
- Keep `annotation.yaml` immutable per cache; store `params.yaml` snapshots.
- Update `vcfcache/public_caches.yaml` when adding/publishing new public aliases (this file is packaged with the wheel).
