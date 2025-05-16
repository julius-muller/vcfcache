# Build stage
FROM python:3.13-slim AS builder

# Install build dependencies and uv
RUN apt update && apt install -y curl && \
    rm -rf /var/lib/apt/lists/* && \
    curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}"

# Only copy what's needed for building
WORKDIR /build
COPY pyproject.toml README.md ./
COPY vcfstash ./vcfstash/

# Build the package (not in development mode)
RUN uv venv /build/venv && \
    . /build/venv/bin/activate && \
    uv pip install --no-cache-dir hatchling && \
    uv pip install --no-cache-dir . && \
    rm -rf ~/.cache

# Runtime stage
FROM python:3.13-slim

# Only install absolutely required runtime dependencies
RUN apt update && apt install -y \
    openjdk-17-jre-headless --no-install-recommends && \
    rm -rf /var/lib/apt/lists/*

# Copy only the built package without development files
WORKDIR /app
COPY --from=builder /build/venv /app/venv
COPY tools /app/tools
COPY resources /app/resources

# Make bcftools executable
RUN chmod +x /app/tools/bcftools

# Set environment variables
ENV PATH="/app/venv/bin:${PATH}"
ENV VCFSTASH_ROOT="/app"

# Create the executable
RUN echo '#!/bin/bash' > /usr/local/bin/vcfstash && \
    echo 'source /app/venv/bin/activate' >> /usr/local/bin/vcfstash && \
    echo 'python -m vcfstash.cli "$@"' >> /usr/local/bin/vcfstash && \
    chmod +x /usr/local/bin/vcfstash

ENTRYPOINT ["/usr/local/bin/vcfstash"]
CMD ["--help"]