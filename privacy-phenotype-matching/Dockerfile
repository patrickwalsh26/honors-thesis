FROM python:3.11-slim

# System deps for cryptography wheel build + unzip for the phenopacket-store
# release archive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        unzip \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python deps first to cache the layer
COPY requirements.txt ./
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Copy the project. The corpus archive lives at data/phenopacket_store/ via
# scripts/download_phenopacket_store.sh and is gitignored; the Makefile
# `reproduce` target fetches it on first run.
COPY . .

ENV PYTHONPATH=/app

# Default: produce all figures and result tables.
CMD ["make", "reproduce"]
