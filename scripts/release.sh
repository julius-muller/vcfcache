#!/usr/bin/env bash
set -euo pipefail

# Release automation script for vcfcache
# Usage: ./scripts/release.sh <version>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Parse arguments
VERSION=""

show_help() {
  cat << EOF
Usage: $0 <version>

Release automation script for vcfcache.
After each step, you'll be prompted to continue, skip, or cancel.

Arguments:
  <version>        Version to release (e.g., 0.3.4)

Options:
  -h, --help       Show this help message

Examples:
  $0 0.3.4

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

version_compare() {
  # Compare two semantic versions (format: X.Y.Z)
  # Returns: 0 if equal, 1 if $1 > $2, 2 if $1 < $2
  local ver1=$1
  local ver2=$2

  if [[ "$ver1" == "$ver2" ]]; then
    return 0
  fi

  local IFS=.
  local i ver1_arr=($ver1) ver2_arr=($ver2)

  # Fill empty positions with zeros
  for ((i=${#ver1_arr[@]}; i<${#ver2_arr[@]}; i++)); do
    ver1_arr[i]=0
  done
  for ((i=${#ver2_arr[@]}; i<${#ver1_arr[@]}; i++)); do
    ver2_arr[i]=0
  done

  for ((i=0; i<${#ver1_arr[@]}; i++)); do
    if ((10#${ver1_arr[i]} > 10#${ver2_arr[i]})); then
      return 1
    elif ((10#${ver1_arr[i]} < 10#${ver2_arr[i]})); then
      return 2
    fi
  done

  return 0
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

log "Starting release process for version $VERSION"
log "Current versions: pyproject.toml=$CURRENT_VERSION_PYPROJECT, __init__.py=$CURRENT_VERSION_INIT"
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

# Get latest version from TestPyPI
log "  → Checking TestPyPI for current version..."
TESTPYPI_LATEST=$(curl -s --max-time 10 "https://test.pypi.org/pypi/vcfcache/json" 2>/dev/null | grep -o '"version":"[^"]*"' | head -1 | cut -d'"' -f4 || echo "")

if [[ -n "$TESTPYPI_LATEST" ]]; then
  set +e
  version_compare "$VERSION" "$TESTPYPI_LATEST"
  COMPARE_RESULT=$?
  set -e
  case $COMPARE_RESULT in
    0)
      log "  ✓ Version v$VERSION already on TestPyPI, skipping"
      ;;
    2)
      log "  ✓ Current TestPyPI version (v$TESTPYPI_LATEST) >= target (v$VERSION), skipping"
      ;;
    1)
      if ask_yn_skip "Upload v$VERSION to TestPyPI? (current: v$TESTPYPI_LATEST)"; then
        python -m twine upload --repository testpypi dist/*
        log "  ✓ Uploaded to TestPyPI"
        log "  → Verify at: https://test.pypi.org/project/vcfcache/$VERSION/"
        echo ""
        read -p "Press Enter once TestPyPI verification is complete..."
      else
        log "  ⊘ Skipped TestPyPI upload"
      fi
      ;;
  esac
else
  # No version on TestPyPI yet
  if ask_yn_skip "Upload v$VERSION to TestPyPI? (current: none)"; then
    python -m twine upload --repository testpypi dist/*
    log "  ✓ Uploaded to TestPyPI"
    log "  → Verify at: https://test.pypi.org/project/vcfcache/$VERSION/"
    echo ""
    read -p "Press Enter once TestPyPI verification is complete..."
  else
    log "  ⊘ Skipped TestPyPI upload"
  fi
fi
echo ""

# Step 4: Upload to PyPI
log "Step 4: Upload to PyPI"

# Get latest version from PyPI
log "  → Checking PyPI for current version..."
PYPI_LATEST=$(curl -s --max-time 10 "https://pypi.org/pypi/vcfcache/json" 2>/dev/null | grep -o '"version":"[^"]*"' | head -1 | cut -d'"' -f4 || echo "")

if [[ -n "$PYPI_LATEST" ]]; then
  set +e
  version_compare "$VERSION" "$PYPI_LATEST"
  COMPARE_RESULT=$?
  set -e
  case $COMPARE_RESULT in
    0)
      log "  ✓ Version v$VERSION already on PyPI, skipping"
      ;;
    2)
      log "  ✓ Current PyPI version (v$PYPI_LATEST) >= target (v$VERSION), skipping"
      ;;
    1)
      if ask_yn_skip "Upload v$VERSION to PyPI? (current: v$PYPI_LATEST)"; then
        python -m twine upload dist/*
        log "  ✓ Uploaded to PyPI"
        log "  → Live at: https://pypi.org/project/vcfcache/$VERSION/"
      else
        log "  ⊘ Skipped PyPI upload"
      fi
      ;;
  esac
else
  # No version on PyPI yet
  if ask_yn_skip "Upload v$VERSION to PyPI? (current: none)"; then
    python -m twine upload dist/*
    log "  ✓ Uploaded to PyPI"
    log "  → Live at: https://pypi.org/project/vcfcache/$VERSION/"
  else
    log "  ⊘ Skipped PyPI upload"
  fi
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
    docker push ghcr.io/julius-muller/vcfcache:latest
    log "  ✓ Pushed Docker images"

    # Verify
    log "  → Verifying pushed image..."
    docker pull ghcr.io/julius-muller/vcfcache:latest --quiet
    docker run --rm ghcr.io/julius-muller/vcfcache:latest demo --smoke-test --quiet
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
  --notes-file CHANGELOG.md \
  dist/*
log "  ✓ GitHub release created"
log "  → View at: https://github.com/julius-muller/vcfcache/releases/tag/v$VERSION"

echo ""
log "================================================"
log "Release v$VERSION complete! 🎉"
log "================================================"
