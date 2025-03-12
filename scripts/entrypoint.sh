# scripts/entrypoint.sh
#!/bin/bash
set -e

# Ensure directories exist
mkdir -p /cache /data

# Check if command is provided
if [ "$1" = "" ]; then
    echo "Usage: vepstash <command> [options]"
    exit 1
fi

# Set VEP environment variables
export VEP_PATH=/opt/vep
export VEP_DATA=/opt/vep/.vep

# Execute command
exec "$@"