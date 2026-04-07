#!/bin/bash
# Usage: bash download_2504_bcf.sh
# Mac から実行 (HPC にインターネット接続なし)
# 2504版 whole-genome BCF を HPC にダウンロード

set -euxo pipefail

DEST_DIR="/Users/ito/mnt/hpc/database/1000genomes/chm13/whole_genome"
S3_BASE="s3://human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1"

echo "Downloading 2504 unrelated samples BCF (~11GB)..."
aws s3 cp "${S3_BASE}/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf.gz" \
    "${DEST_DIR}/" --no-sign-request

aws s3 cp "${S3_BASE}/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf.gz.csi" \
    "${DEST_DIR}/" --no-sign-request

echo "Done."
ls -lh "${DEST_DIR}/"
