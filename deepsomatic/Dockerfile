# DeepSomatic Docker Image
# Based on official Google DeepSomatic image
# https://github.com/google/deepsomatic

FROM google/deepsomatic:1.9.0

LABEL maintainer="itoyu8"
LABEL description="DeepSomatic - Deep learning based somatic variant caller for tumor-normal paired sequencing"

# The official image includes:
# - DeepSomatic binary
# - Pre-trained models for WGS, WES, PACBIO, ONT, FFPE
# - Tumor-only models available

# Supported model_type:
#   - Tumor-normal: WGS, WES, PACBIO, ONT, FFPE_WGS, FFPE_WES
#   - Tumor-only: WGS_TUMOR_ONLY, PACBIO_TUMOR_ONLY, ONT_TUMOR_ONLY, FFPE_WGS_TUMOR_ONLY, FFPE_WES_TUMOR_ONLY

WORKDIR /data

ENTRYPOINT ["run_deepsomatic"]
