# Multi-stage Dockerfile for DeepVariant and WhatsHap
# Copyright 2024 - Based on Google DeepVariant Dockerfile
# Combines variant calling (DeepVariant) with phasing (WhatsHap)

ARG FROM_IMAGE=ubuntu:22.04
ARG DV_GPU_BUILD=0
ARG VERSION=1.9.0
ARG PYTHON_VERSION=3.10

# Stage 1: Setup conda environment
FROM continuumio/miniconda3:23.3.1-0 AS conda_setup
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set offline false && \
    conda install -y bcftools=1.15 samtools=1.15 whatshap && \
    conda clean -a

# Stage 2: Build DeepVariant binaries
FROM ${FROM_IMAGE} AS builder
ARG DV_GPU_BUILD
ARG VERSION
ARG PYTHON_VERSION

LABEL maintainer="itoyu8"
LABEL description="Integrated container for DeepVariant and WhatsHap"

# Copy conda environment
COPY --from=conda_setup /opt/conda /opt/conda
ENV PATH="/opt/conda/bin:${PATH}"

# Install build dependencies
RUN apt-get -qq update && \
    apt-get -qq install -y \
    build-essential \
    curl \
    git \
    python3-dev \
    python3-pip \
    wget \
    zip \
    && apt-get clean

# Install Bazel (required for building DeepVariant/DeepSomatic)
RUN wget https://github.com/bazelbuild/bazel/releases/download/6.1.1/bazel-6.1.1-installer-linux-x86_64.sh && \
    bash bazel-6.1.1-installer-linux-x86_64.sh && \
    rm bazel-6.1.1-installer-linux-x86_64.sh

# Clone and build DeepVariant
WORKDIR /opt
RUN git clone https://github.com/google/deepvariant.git && \
    cd deepvariant && \
    git checkout r${VERSION}

WORKDIR /opt/deepvariant
RUN ./build-prereq.sh && \
    DV_GPU_BUILD=${DV_GPU_BUILD} ./build_release_binaries.sh

# Stage 3: Runtime image
FROM ${FROM_IMAGE}
ARG DV_GPU_BUILD
ARG VERSION
ARG PYTHON_VERSION

# Copy conda environment
COPY --from=conda_setup /opt/conda /opt/conda
ENV PATH="/opt/conda/bin:${PATH}"

# Copy built binaries from builder
COPY --from=builder /opt/deepvariant/bazel-bin/deepvariant /opt/deepvariant/bin/
COPY --from=builder /opt/deepvariant/scripts /opt/deepvariant/scripts/
COPY --from=builder /opt/deepvariant/settings.sh /opt/deepvariant/

# Install runtime dependencies
RUN apt-get -qq update && \
    apt-get -qq install -y \
    parallel \
    python${PYTHON_VERSION} \
    python3-pip \
    wget \
    curl \
    tabix \
    && apt-get clean

# Set Python 3.10 as default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python${PYTHON_VERSION} 1 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python${PYTHON_VERSION} 1

# Install Python dependencies for DeepVariant/DeepSomatic
RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install \
    absl-py \
    contextlib2 \
    numpy \
    pillow \
    tensorflow \
    tensorflow-io-gcs-filesystem

# Download DeepVariant models
WORKDIR /opt/models/deepvariant
RUN curl https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wgs_standard/model.ckpt.data-00000-of-00001 > wgs.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wgs_standard/model.ckpt.index > wgs.model.ckpt.index && \
    curl https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wes_standard/model.ckpt.data-00000-of-00001 > wes.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wes_standard/model.ckpt.index > wes.model.ckpt.index && \
    curl https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-pacbio_standard/model.ckpt.data-00000-of-00001 > pacbio.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-pacbio_standard/model.ckpt.index > pacbio.model.ckpt.index && \
    curl https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-ont_r104_e81_sup_g5015/model.ckpt.data-00000-of-00001 > ont.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-ont_r104_e81_sup_g5015/model.ckpt.index > ont.model.ckpt.index

# Create shell wrappers for DeepVariant
WORKDIR /opt/deepvariant/bin
RUN for bin in make_examples call_variants postprocess_variants; do \
      echo '#!/bin/bash' > ${bin} && \
      echo 'set -euo pipefail' >> ${bin} && \
      echo 'python3 /opt/deepvariant/bin/'${bin}'.zip "$@"' >> ${bin} && \
      chmod +x ${bin}; \
    done

RUN echo '#!/bin/bash' > run_deepvariant && \
    echo 'set -euo pipefail' >> run_deepvariant && \
    echo 'python3 /opt/deepvariant/scripts/run_deepvariant.py "$@"' >> run_deepvariant && \
    chmod +x run_deepvariant

# Add binaries to PATH
ENV PATH="/opt/deepvariant/bin:${PATH}"

# Verify installations
RUN python3 --version && \
    bcftools --version && \
    samtools --version && \
    whatshap --version && \
    tabix --version

WORKDIR /data

# Default command shows help
CMD ["/bin/bash", "-c", "echo 'DeepVariant and WhatsHap container'; echo ''; echo 'Available commands:'; echo '  - run_deepvariant'; echo '  - whatshap'; echo '  - bcftools'; echo '  - samtools'; echo '  - tabix'"]
