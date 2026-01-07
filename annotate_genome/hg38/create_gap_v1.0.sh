#!/bin/bash
# Create hg38_gap_v1.0.bed from UCSC centromeres and gap annotations
# Usage: ./create_gap_bed.sh
#
# Source: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/
# - centromeres.txt: Centromere model sequence positions
# - gap.txt: Assembly gap annotations

SCRIPT_DIR=$(dirname "$0")
INTERMEDIATE="${SCRIPT_DIR}/intermediate"
CENTROMERES="${INTERMEDIATE}/centromeres.txt"
GAP="${INTERMEDIATE}/gap.txt"
OUTPUT="${SCRIPT_DIR}/hg38_gap_v1.0.bed"

# 1. Process centromeres: min-max per chromosome
# centromeres.txt format: bin chr start end name
awk -F'\t' 'BEGIN {OFS="\t"}
{
    chr = $2
    start = $3
    end = $4
    if (!(chr in min) || start < min[chr]) min[chr] = start
    if (!(chr in max) || end > max[chr]) max[chr] = end
}
END {
    for (chr in min) print chr, min[chr], max[chr], "centromere"
}' "${CENTROMERES}" | sort -k1,1V -k2,2n > /tmp/hg38_centromeres.bed

# 2. Extract gap regions: telomere, heterochromatin, short_arm
# gap.txt format: bin chr start end ix N size type bridge
awk -F'\t' 'BEGIN {OFS="\t"} $8=="telomere" || $8=="heterochromatin" || $8=="short_arm" {
    print $2, $3, $4, $8
}' "${GAP}" | sort -k1,1V -k2,2n > /tmp/hg38_gaps.bed

# 3. Combine, sort, and format as BED9
{
    cat /tmp/hg38_centromeres.bed
    cat /tmp/hg38_gaps.bed
} | sort -k1,1V -k2,2n | \
    awk -F'\t' 'BEGIN {OFS="\t"}
    {
        # Assign color based on region type
        if ($4 == "centromere") color = "224,0,0"
        else if ($4 == "telomere") color = "0,0,224"
        else if ($4 == "heterochromatin") color = "128,128,128"
        else if ($4 == "short_arm") color = "224,224,0"
        else color = "0,0,0"
        print $1, $2, $3, $4, 100, ".", $2, $3, color
    }' > "${OUTPUT}.tmp"

# 4. Add header
echo 'track name="hg38_exclude" description="hg38 centromere/telomere/heterochromatin/short_arm regions for VCF filtering" itemRgb="On" visibility="1"' | \
    cat - "${OUTPUT}.tmp" > "${OUTPUT}"

rm -f "${OUTPUT}.tmp" /tmp/hg38_centromeres.bed /tmp/hg38_gaps.bed

echo "Created: ${OUTPUT}"
wc -l "${OUTPUT}"
