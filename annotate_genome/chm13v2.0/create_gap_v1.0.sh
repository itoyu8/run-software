#!/bin/bash
# Create chm13v2.0_gap_v1.0.bed from censat and telomere annotations
# Usage: ./create_gap_bed.sh

SCRIPT_DIR=$(dirname "$0")
CENSAT="${SCRIPT_DIR}/chm13v2.0_censat_v2.0.bed"
TELOMERE="${SCRIPT_DIR}/chm13v2.0_telomere.bed"
OUTPUT="${SCRIPT_DIR}/chm13v2.0_gap_v1.0.bed"

# Extract hor/hsat/bsat/gsat from censat, add telomere, sort, and merge
{
    # Extract satellite regions (hor, hsat, bsat, gsat) from censat
    grep -E 'hor|hsat|bsat|gsat' "${CENSAT}" | cut -f1-3

    # Add telomere regions
    cat "${TELOMERE}"
} | \
    sort -k1,1V -k2,2n | \
    awk -F'\t' '
    BEGIN { OFS="\t" }
    NR==1 { chr=$1; start=$2; end=$3; next }
    {
        if ($1 == chr && $2 <= end) {
            # Overlapping or adjacent -> merge
            if ($3 > end) end = $3
        } else {
            # New region
            print chr, start, end
            chr = $1; start = $2; end = $3
        }
    }
    END { print chr, start, end }
    ' | \
    awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3, "censat", 100, ".", $2, $3, "224,0,0"}' > "${OUTPUT}.tmp"

# Add header
echo 'track name="chm13_exclude" description="chm13 centromeric/satellite/telomere regions for VCF filtering" itemRgb="On" visibility="1"' | \
    cat - "${OUTPUT}.tmp" > "${OUTPUT}"

rm -f "${OUTPUT}.tmp"

echo "Created: ${OUTPUT}"
wc -l "${OUTPUT}"
