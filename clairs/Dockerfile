# ClairS Docker Image
# Based on official HKU-BAL ClairS image
# https://github.com/HKU-BAL/ClairS

FROM hkubal/clairs:latest

LABEL maintainer="itoyu8"
LABEL description="ClairS - Somatic small variant caller for tumor-normal paired long-read sequencing"

# The official image includes:
# - Python 3.9+
# - TensorFlow
# - Samtools
# - WhatsHap
# - Pre-trained models for ONT and PacBio

# Models are located at /opt/models/ in the container
# Available platforms:
#   - ONT: ont_r10_dorado_sup_5khz, ont_r10_dorado_hac_5khz, etc.
#   - PacBio HiFi: hifi_revio, hifi_sequel2
#   - Illumina: ilmn

WORKDIR /data

ENTRYPOINT ["/opt/bin/run_clairs"]
