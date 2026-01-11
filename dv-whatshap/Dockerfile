# DeepVariant + WhatsHap Docker Image
# Based on official Google DeepVariant image
# https://github.com/google/deepvariant

FROM google/deepvariant:1.9.0

LABEL maintainer="itoyu8"
LABEL description="DeepVariant with WhatsHap for variant calling and phasing"

# The official image includes:
# - DeepVariant binary (run_deepvariant)
# - Pre-trained models for WGS, WES, PACBIO, ONT_R104

# Add WhatsHap for phasing
RUN apt-get update && \
    apt-get install -y python3-pip && \
    pip3 install whatshap && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /data
