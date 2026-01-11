FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libncurses5-dev \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# BWA 0.7.19 (from source archive)
RUN wget https://github.com/lh3/bwa/archive/refs/tags/v0.7.19.tar.gz -O bwa-0.7.19.tar.gz && \
    tar -xzf bwa-0.7.19.tar.gz && \
    cd bwa-0.7.19 && \
    make && \
    cp bwa /usr/local/bin/

# Samtools 1.19 (includes htslib)
RUN wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar -xjf samtools-1.19.tar.bz2 && \
    cd samtools-1.19 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install

# BCFtools 1.19
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 && \
    tar -xjf bcftools-1.19.tar.bz2 && \
    cd bcftools-1.19 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install

# Minimap2 2.28
RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2 && \
    tar -xjf minimap2-2.28.tar.bz2 && \
    cd minimap2-2.28 && \
    make && \
    cp minimap2 /usr/local/bin/

# gfatools (latest)
RUN wget https://github.com/lh3/gfatools/archive/refs/heads/master.tar.gz -O gfatools.tar.gz && \
    tar -xzf gfatools.tar.gz && \
    cd gfatools-master && \
    make && \
    cp gfatools /usr/local/bin/

# seqtk (lh3, latest)
RUN wget https://github.com/lh3/seqtk/archive/refs/heads/master.tar.gz -O seqtk.tar.gz && \
    tar -xzf seqtk.tar.gz && \
    cd seqtk-master && \
    make && \
    cp seqtk /usr/local/bin/

# seqkit (pre-built binary, Go-based)
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz && \
    tar -xzf seqkit_linux_amd64.tar.gz && \
    cp seqkit /usr/local/bin/

# Runtime stage
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    libncurses6 \
    && rm -rf /var/lib/apt/lists/*

# Copy binaries from builder
COPY --from=builder /usr/local/bin/bwa /usr/local/bin/
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/
COPY --from=builder /usr/local/bin/bcftools /usr/local/bin/
COPY --from=builder /usr/local/bin/minimap2 /usr/local/bin/
COPY --from=builder /usr/local/bin/gfatools /usr/local/bin/
COPY --from=builder /usr/local/bin/seqtk /usr/local/bin/
COPY --from=builder /usr/local/bin/seqkit /usr/local/bin/

# Verify installations
RUN bwa 2>&1 | head -3 && \
    samtools --version | head -2 && \
    bcftools --version | head -2 && \
    minimap2 --version && \
    gfatools version && \
    seqtk 2>&1 | head -3 && \
    seqkit version

WORKDIR /data

CMD ["/bin/bash"]
