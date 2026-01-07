# Stage 1: Build environment
FROM rust:1.75-slim-bookworm AS builder

WORKDIR /build

# Clone and build rukki
RUN apt-get update && apt-get install -y git && \
    git clone https://github.com/marbl/rukki.git . && \
    cargo build --release

# Stage 2: Runtime environment
FROM debian:bookworm-slim

# Install Python and necessary libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgcc-s1 \
    python3 \
    python3-pip \
    && pip3 install --no-cache-dir --break-system-packages networkx numpy \
    && rm -rf /var/lib/apt/lists/*

# Copy binary from builder
COPY --from=builder /build/target/release/rukki /usr/local/bin/

# Copy verkko scripts for paths2fasta conversion
COPY verkko_scripts_v2.2.1 /opt/verkko_scripts/

ENTRYPOINT ["rukki"]
