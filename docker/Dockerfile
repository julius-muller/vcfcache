FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    python3 python3-pip openjdk-17-jre-headless procps curl \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Create Python symlink
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install uv package manager and add it to PATH
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}"

# Set working directory
WORKDIR /app

# Copy application files
COPY . /app

# Make bcftools executable
RUN chmod +x /app/tools/bcftools

# Create venv and install dependencies
RUN /root/.local/bin/uv venv /app/.venv && \
    . /app/.venv/bin/activate && \
    /root/.local/bin/uv pip install hatchling editables && \
    /root/.local/bin/uv pip install -e ".[parquet]"

# Set environment variables for runtime
ENV PATH="/app/.venv/bin:${PATH}"
ENV VCFSTASH_ROOT="/app"

# Use ENTRYPOINT with the proper command to run the CLI
# Looking at cli.py, we should run the main() function directly
ENTRYPOINT ["python", "-c", "from vcfstash.cli import main; main()"]

# Default command is empty since arguments will be passed to the entrypoint
CMD []