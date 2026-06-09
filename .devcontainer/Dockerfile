FROM mcr.microsoft.com/devcontainers/python:3.11-bookworm

# System dependencies: wget (downloads) and clustalo (alignment)
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        wget \
        clustalo && \
    rm -rf /var/lib/apt/lists/*

# Install MMseqs2 (SSE4.1 build) and expose it on the PATH as `mmseqs`
RUN wget -q https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz -O /tmp/mmseqs.tar.gz && \
    tar -xzf /tmp/mmseqs.tar.gz -C /tmp && \
    cp /tmp/mmseqs/bin/mmseqs /usr/local/bin/mmseqs && \
    chmod +x /usr/local/bin/mmseqs && \
    rm -rf /tmp/mmseqs /tmp/mmseqs.tar.gz

# uv for dependency management
RUN pip install uv

USER vscode
