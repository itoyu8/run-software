# Dockerfile for rasusa - Rust-based FASTQ downsampling tool
FROM rust:1.85-slim AS builder

# Install dependencies for building
RUN apt-get update && apt-get install -y \
    curl \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install rasusa using cargo
RUN cargo install rasusa

# Final stage - minimal image
FROM debian:bookworm-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy rasusa binary from builder
COPY --from=builder /usr/local/cargo/bin/rasusa /usr/local/bin/rasusa

# Set working directory
WORKDIR /data

# Default command - show help
CMD ["rasusa", "--help"]
