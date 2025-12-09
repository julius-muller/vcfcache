# Releasing `vcfcache`

This document describes the steps to create a new `vcfcache` release and publish
it to PyPI. The flow is:

1. Bump version
2. Update changelog & docs
3. Build & upload to **TestPyPI**
4. Smoke-test from TestPyPI
5. Upload to **PyPI**
6. Tag and create a GitHub release

> Notes:
> - The project is built with `hatchling` via `pyproject.toml`.
> - Distributions are built with `python -m build`.
> - Uploads are done with `twine`.
> - Python version and metadata live in `pyproject.toml`. :contentReference[oaicite:0]{index=0}

---

## 0. Prerequisites

- Python >= 3.11 installed
- `uv` or `pip` to manage the local environment
- Accounts on:
  - https://pypi.org
  - https://test.pypi.org
- API tokens created on both:
  - Go to “Account settings → API tokens”
  - Scope: whole account is fine for this project
- Locally:
  ```bash
  uv venv .venv
  source .venv/bin/activate
  uv pip install build twine
  # or: pip install build twine
