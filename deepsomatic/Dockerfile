# Multi-stage Dockerfile for DeepSomatic
# Copyright 2024 - Based on Google DeepSomatic Dockerfile
# Somatic variant calling for tumor/normal samples

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
    conda install -y bcftools=1.15 samtools=1.15 && \
    conda clean -a

# Stage 2: Build DeepSomatic binaries
FROM ${FROM_IMAGE} AS builder
ARG DV_GPU_BUILD
ARG VERSION
ARG PYTHON_VERSION

LABEL maintainer="itoyu8"
LABEL description="Container for DeepSomatic somatic variant calling"

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

# Install Bazel (required for building DeepSomatic)
RUN wget https://github.com/bazelbuild/bazel/releases/download/6.1.1/bazel-6.1.1-installer-linux-x86_64.sh && \
    bash bazel-6.1.1-installer-linux-x86_64.sh && \
    rm bazel-6.1.1-installer-linux-x86_64.sh

# Clone and build DeepVariant (DeepSomatic is part of it)
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
COPY --from=builder /opt/deepvariant/bazel-bin/deepsomatic /opt/deepsomatic/bin/
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

# Install Python dependencies for DeepSomatic
RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install \
    absl-py \
    contextlib2 \
    numpy \
    pillow \
    tensorflow \
    tensorflow-io-gcs-filesystem

# Download DeepSomatic models
WORKDIR /opt/models/deepsomatic

# WGS models (somatic and tumor-only)
RUN curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-wgs_somatic/model.ckpt.data-00000-of-00001 > wgs_somatic.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-wgs_somatic/model.ckpt.index > wgs_somatic.model.ckpt.index && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-wgs_tumor_only/model.ckpt.data-00000-of-00001 > wgs_tumor_only.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-wgs_tumor_only/model.ckpt.index > wgs_tumor_only.model.ckpt.index

# WES models (somatic and tumor-only)
RUN curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-wes_somatic/model.ckpt.data-00000-of-00001 > wes_somatic.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-wes_somatic/model.ckpt.index > wes_somatic.model.ckpt.index && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-wes_tumor_only/model.ckpt.data-00000-of-00001 > wes_tumor_only.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-wes_tumor_only/model.ckpt.index > wes_tumor_only.model.ckpt.index

# PacBio models (somatic and tumor-only)
RUN curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-pacbio_somatic/model.ckpt.data-00000-of-00001 > pacbio_somatic.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-pacbio_somatic/model.ckpt.index > pacbio_somatic.model.ckpt.index && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-pacbio_tumor_only/model.ckpt.data-00000-of-00001 > pacbio_tumor_only.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-pacbio_tumor_only/model.ckpt.index > pacbio_tumor_only.model.ckpt.index

# ONT models (somatic and tumor-only)
RUN curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-ont_r104_e81_sup_g5015_somatic/model.ckpt.data-00000-of-00001 > ont_somatic.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-ont_r104_e81_sup_g5015_somatic/model.ckpt.index > ont_somatic.model.ckpt.index && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-ont_r104_e81_sup_g5015_tumor_only/model.ckpt.data-00000-of-00001 > ont_tumor_only.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-ont_r104_e81_sup_g5015_tumor_only/model.ckpt.index > ont_tumor_only.model.ckpt.index

# FFPE models (WGS and WES)
RUN curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-ffpe_wgs_somatic/model.ckpt.data-00000-of-00001 > ffpe_wgs_somatic.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-ffpe_wgs_somatic/model.ckpt.index > ffpe_wgs_somatic.model.ckpt.index && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-ffpe_wes_somatic/model.ckpt.data-00000-of-00001 > ffpe_wes_somatic.model.ckpt.data-00000-of-00001 && \
    curl https://storage.googleapis.com/deepvariant/models/DeepSomatic/${VERSION}/DeepSomatic-inception_v3-${VERSION}+data-ffpe_wes_somatic/model.ckpt.index > ffpe_wes_somatic.model.ckpt.index

# Create shell wrappers for DeepSomatic
WORKDIR /opt/deepsomatic/bin
RUN for bin in make_examples_somatic call_variants postprocess_variants; do \
      echo '#!/bin/bash' > ${bin} && \
      echo 'set -euo pipefail' >> ${bin} && \
      echo 'python3 /opt/deepsomatic/bin/'${bin}'.zip "$@"' >> ${bin} && \
      chmod +x ${bin}; \
    done

RUN echo '#!/bin/bash' > run_deepsomatic && \
    echo 'set -euo pipefail' >> run_deepsomatic && \
    echo 'python3 /opt/deepvariant/scripts/run_deepsomatic.py "$@"' >> run_deepsomatic && \
    chmod +x run_deepsomatic

# Add binaries to PATH
ENV PATH="/opt/deepsomatic/bin:${PATH}"

# Verify installations
RUN python3 --version && \
    bcftools --version && \
    samtools --version && \
    tabix --version

WORKDIR /data

# Default command shows help
CMD ["/bin/bash", "-c", "echo 'DeepSomatic container'; echo ''; echo 'Available commands:'; echo '  - run_deepsomatic'; echo '  - make_examples_somatic'; echo '  - call_variants'; echo '  - postprocess_variants'; echo '  - bcftools'; echo '  - samtools'; echo '  - tabix'"]
