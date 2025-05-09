#!/bin/bash
set -e

# Default directories
DEFAULT_REFERENCE_DIR="./reference"
DEFAULT_DATA_DIR="./data"
DEFAULT_CACHE_DIR="./cache"

# Use environment variables if set, otherwise use defaults
REFERENCE_DIR=${REFERENCE_DIR:-$DEFAULT_REFERENCE_DIR}
DATA_DIR=${DATA_DIR:-$DEFAULT_DATA_DIR}
CACHE_DIR=${CACHE_DIR:-$DEFAULT_CACHE_DIR}

# Export variables for docker-compose
export REFERENCE_DIR
export DATA_DIR
export CACHE_DIR

# Run docker with the provided arguments
docker run --rm \
  -v ${REFERENCE_DIR}:/reference:ro \
  -v ${DATA_DIR}:/data \
  -v ${CACHE_DIR}:/cache \
  -e REFERENCE_PATH=/reference \
  vcfstash:latest "$@"