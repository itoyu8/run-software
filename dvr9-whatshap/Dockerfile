# PEPPER-Margin-DeepVariant + WhatsHap Docker Image
# For ONT R9.4.1 data variant calling and phasing
# https://github.com/kishwarshafin/pepper

FROM kishwars/pepper_deepvariant:r0.8

LABEL maintainer="itoyu8"
LABEL description="PEPPER-Margin-DeepVariant with WhatsHap for ONT R9.4.1 variant calling and phasing"

# The official image includes:
# - PEPPER-Margin-DeepVariant pipeline (run_pepper_margin_deepvariant)
# - Pre-trained models for ONT R9.4.1 (--ont_r9_guppy5_sup)

# Add WhatsHap for additional phasing capabilities
RUN apt-get update && \
    apt-get install -y python3-pip && \
    pip3 install whatshap && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /data
