#!/usr/bin/env bash
set -euo pipefail

# Release automation script for vcfcache
# Usage: ./scripts/release.sh <version>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Parse arguments
VERSION=""
FORCE_GH_PRERELEASE=0

show_help() {
  cat << EOF
Usage: $0 <version>

Release automation script for vcfcache.
After each step, you'll be prompted to continue, skip, or cancel.

Arguments:
  <version>        Version to release (e.g., 0.4.0, 0.4.0b0, 0.4.0rc1)

Options:
  -h, --help       Show this help message
  --github-prerelease  Mark the GitHub Release as pre-release (even if version is not a PEP 440 pre-release)

Examples:
  $0 0.3.4
  $0 0.4.0b0
  $0 0.4.0 --github-prerelease

Interactive prompts:
  y - Yes, proceed with this step
  n - No, cancel the release and exit
  s - Skip this step and continue
EOF
}

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      show_help
      exit 0
      ;;
    --github-prerelease)
      FORCE_GH_PRERELEASE=1
      shift
      ;;
    -*)
      echo "Unknown option: $1"
      show_help
      exit 1
      ;;
    *)
      if [[ -z "$VERSION" ]]; then
        VERSION="$1"
      else
        echo "Error: Unexpected argument: $1"
        show_help
        exit 1
      fi
      shift
      ;;
  esac
done

if [[ -z "$VERSION" ]]; then
  echo "Error: Version is required"
  show_help
  exit 1
fi

log() { echo "[$(date +%H:%M:%S)] $*"; }

is_pep440_prerelease() {
  # Prefer a real PEP 440 check (packaging); fall back to a heuristic.
  local v="$1"
  if python -c "from packaging.version import Version; import sys; print('1' if Version(sys.argv[1]).is_prerelease else '0')" "$v" >/dev/null 2>&1; then
    python -c "from packaging.version import Version; import sys; print('1' if Version(sys.argv[1]).is_prerelease else '0')" "$v"
    return 0
  fi
  if [[ "$v" =~ (a|b|rc)[0-9]+$ ]]; then
    echo "1"
  else
    echo "0"
  fi
}

ask_yn_skip() {
  local prompt="$1"
  while true; do
    read -p "$prompt (y=yes, n=cancel, s=skip) " -n 1 -r
    echo
    case $REPLY in
      [Yy])
        return 0  # yes
        ;;
      [Nn])
        log "Cancelled by user"
        exit 1
        ;;
      [Ss])
        return 1  # skip
        ;;
      *)
        echo "Invalid input. Please enter y (yes), n (cancel), or s (skip)"
        ;;
    esac
  done
}

cd "$PROJECT_ROOT"

# Get current versions from files
CURRENT_VERSION_PYPROJECT=$(grep '^version = ' pyproject.toml | sed 's/version = "\(.*\)"/\1/')
CURRENT_VERSION_INIT=$(grep '^__version__ = ' vcfcache/__init__.py | sed 's/__version__ = "\(.*\)"/\1/')
IS_PRERELEASE="$(is_pep440_prerelease "$VERSION")"
GH_PRERELEASE_FLAG=""
if [[ "$IS_PRERELEASE" == "1" ]] || [[ "$FORCE_GH_PRERELEASE" == "1" ]]; then
  GH_PRERELEASE_FLAG="--prerelease"
fi

log "Starting release process for version $VERSION"
log "Current versions: pyproject.toml=$CURRENT_VERSION_PYPROJECT, __init__.py=$CURRENT_VERSION_INIT"
if [[ "$IS_PRERELEASE" == "1" ]]; then
  log "Detected pre-release version (PEP 440): $VERSION"
fi
if [[ "$FORCE_GH_PRERELEASE" == "1" ]] && [[ "$IS_PRERELEASE" != "1" ]]; then
  log "GitHub Release will be marked as pre-release (--github-prerelease)"
fi
echo ""

# Step 1: Update version in files
log "Step 1: Update version in pyproject.toml, __init__.py, CHANGELOG.md"

# Check if version is already updated
if [[ "$CURRENT_VERSION_PYPROJECT" == "$VERSION" ]] && [[ "$CURRENT_VERSION_INIT" == "$VERSION" ]]; then
  log "  ✓ Version already set to $VERSION, skipping"
else
  if ask_yn_skip "Proceed with version update? ($CURRENT_VERSION_PYPROJECT -> $VERSION)"; then

    # Update pyproject.toml
    if [[ -f "pyproject.toml" ]]; then
      sed -i "s/^version = .*/version = \"$VERSION\"/" pyproject.toml
      log "  ✓ Updated pyproject.toml"
    else
      log "  ✗ pyproject.toml not found"
      exit 1
    fi

    # Update __init__.py
    if [[ -f "vcfcache/__init__.py" ]]; then
      sed -i "s/^__version__ = .*/__version__ = \"$VERSION\"/" vcfcache/__init__.py
      log "  ✓ Updated vcfcache/__init__.py"
    else
      log "  ✗ vcfcache/__init__.py not found"
      exit 1
    fi

    # Reminder to update CHANGELOG.md manually
    log "  ⚠ Please manually update CHANGELOG.md with:"
    log "    ## [$VERSION] - $(date +%Y-%m-%d)"
    log "    ### Added/Changed/Fixed"
    log "    - ..."
    echo ""
    read -p "Press Enter once CHANGELOG.md is updated..."
    log "  ✓ Version update complete"
  else
    log "  ⊘ Skipped version update"
  fi
fi
echo ""

# Step 2: Build and test locally
log "Step 2: Build and test locally"

# Check if package is already built
if [[ -f "dist/vcfcache-${VERSION}-py3-none-any.whl" ]] && [[ -f "dist/vcfcache-${VERSION}.tar.gz" ]]; then
  log "  ✓ Package v$VERSION already built, skipping"
else
  if ask_yn_skip "Build package v$VERSION and run tests?"; then
    rm -rf dist/ build/ *.egg-info
    python -m build
    log "  ✓ Built package"

    log "  → Installing in temporary venv and running tests..."
    uv venv /tmp/vcfcache-release-test
    source /tmp/vcfcache-release-test/bin/activate
    uv pip install "dist/vcfcache-${VERSION}-py3-none-any.whl[dev]"

    log "  → Running smoke test..."
    vcfcache demo --smoke-test --quiet

    log "  → Running test suite..."
    python -m pytest tests -q

    deactivate
    rm -rf /tmp/vcfcache-release-test

    log "  ✓ Local tests passed"
  else
    log "  ⊘ Skipped build and test"
  fi
fi
echo ""

# Step 3: Upload to TestPyPI
log "Step 3: Upload to TestPyPI"

if [[ "$IS_PRERELEASE" == "1" ]]; then
  log "  → Pre-release detected: TestPyPI is recommended for betas/RCs"
fi
if ask_yn_skip "Upload v$VERSION to TestPyPI?"; then
  python -m twine upload --repository testpypi dist/*
  log "  ✓ Uploaded to TestPyPI"
  log "  → Verify at: https://test.pypi.org/project/vcfcache/$VERSION/"
  echo ""
  read -p "Press Enter once TestPyPI verification is complete..."
else
  log "  ⊘ Skipped TestPyPI upload"
fi
echo ""

# Step 4: Upload to PyPI
log "Step 4: Upload to PyPI"

if [[ "$IS_PRERELEASE" == "1" ]]; then
  log "  ⚠ Pre-release detected: skipping PyPI upload is usually recommended"
fi
if ask_yn_skip "Upload v$VERSION to PyPI?"; then
  python -m twine upload dist/*
  log "  ✓ Uploaded to PyPI"
  log "  → Live at: https://pypi.org/project/vcfcache/$VERSION/"
else
  log "  ⊘ Skipped PyPI upload"
fi
echo ""

# Step 5: Docker build and push
log "Step 5: Build and push Docker image"

# Ask user if they want to build/push Docker image
if ! ask_yn_skip "Build Docker image for v$VERSION?"; then
  log "  ⊘ Skipped Docker build"
else
  # Get latest Docker image version from GHCR (optional, for informational purposes)
  log "  → Checking GHCR for current version..."
  DOCKER_LATEST=""
  if docker pull ghcr.io/julius-muller/vcfcache:latest --quiet 2>/dev/null; then
    DOCKER_LATEST=$(docker run --rm ghcr.io/julius-muller/vcfcache:latest --version 2>/dev/null | head -1 || echo "")
    if [[ -n "$DOCKER_LATEST" ]]; then
      log "  → Current GHCR version: v$DOCKER_LATEST"
    fi
  fi

  # Build image
  log "  → Building Docker image..."
  ./scripts/local-build/build-and-push-final.sh --skip-push --force

  # Tag with version
  log "  → Tagging as v$VERSION..."
  docker tag ghcr.io/julius-muller/vcfcache:latest "ghcr.io/julius-muller/vcfcache:v$VERSION"

  # Push
  if ask_yn_skip "Push Docker images to GHCR?"; then
    docker push "ghcr.io/julius-muller/vcfcache:v$VERSION"
    if [[ "$IS_PRERELEASE" != "1" ]]; then
      docker push ghcr.io/julius-muller/vcfcache:latest
    else
      log "  → Pre-release detected: not pushing :latest"
    fi
    log "  ✓ Pushed Docker images"

    # Verify
    log "  → Verifying pushed image..."
    docker pull "ghcr.io/julius-muller/vcfcache:v$VERSION" --quiet
    docker run --rm "ghcr.io/julius-muller/vcfcache:v$VERSION" demo --smoke-test --quiet
    log "  ✓ Docker image verified"
  else
    log "  ⊘ Skipped Docker push (images tagged locally)"
  fi
fi
echo ""

# Step 6: GitHub release
log "Step 6: Create GitHub release"

# Check if release already exists
if gh release view "v$VERSION" >/dev/null 2>&1; then
  log "  ✓ GitHub release v$VERSION already exists, skipping"
  exit 0
fi

# Ask user if they want to create the release
if ! ask_yn_skip "Create GitHub release v$VERSION with artifacts?"; then
  log "  ⊘ Skipped GitHub release creation"
  exit 0
fi

# Check if git tag exists
if git rev-parse "v$VERSION" >/dev/null 2>&1; then
  log "  ✓ Git tag v$VERSION already exists"
else
  # Create and push git tag
  log "  → Creating git tag v$VERSION..."
  git tag -a "v$VERSION" -m "Release v$VERSION"
  git push origin "v$VERSION"
  log "  ✓ Created and pushed git tag v$VERSION"
fi

# Create GitHub release
log "  → Creating GitHub release..."
gh release create "v$VERSION" \
  --title "vcfcache v$VERSION" \
  $GH_PRERELEASE_FLAG \
  --notes-file CHANGELOG.md \
  dist/*
log "  ✓ GitHub release created"
log "  → View at: https://github.com/julius-muller/vcfcache/releases/tag/v$VERSION"

echo ""
log "================================================"
log "Release v$VERSION complete! 🎉"
log "================================================"
