#!/bin/bash
# Create hg38_gap_v1.1.bed
# Combines: centromere/heterochromatin/short_arm/telomere (from gap.txt), segdup (>10kb, >99% identity)
# Output: BED9 format with track header, categories distinguished by name and color
# Usage: ./create_exclude_bed.sh

set -euxo pipefail

SCRIPT_DIR=$(dirname "$(realpath "$0")")
cd "$SCRIPT_DIR"

# Source files
GAP_BED="hg38_gap_v1.0.bed"
SEGDUP_RAW="../../data/genome_stratifications/ucsc_segdup/hg38/genomicSuperDups.txt.gz"
OUTPUT="hg38_gap_v1.1.bed"

# Colors (RGB) - same as v1.0 plus segdup
COLOR_CENTROMERE="224,0,0"      # Red
COLOR_TELOMERE="0,0,224"        # Blue
COLOR_HETEROCHROMATIN="128,128,128"  # Gray
COLOR_SHORT_ARM="224,224,0"     # Yellow
COLOR_SEGDUP="0,128,0"          # Green

# 1. Extract existing gap regions (already in BED9 format with categories)
awk 'NR>1 {print}' "$GAP_BED" > gap_regions.tmp

# 2. Segdup: filter by length>=10kb, identity>=99%, merge, output as BED9
gzip -dc "$SEGDUP_RAW" | awk -F'\t' '
BEGIN {OFS="\t"}
{
    chrom = $2; start = $3; end = $4
    len = end - start
    fracMatch = $27
    if (len >= 10000 && fracMatch >= 0.99) {
        print chrom, start, end
    }
}' | sort -k1,1V -k2,2n | awk '
BEGIN {OFS="\t"}
NR==1 {chr=$1; start=$2; end=$3; next}
{
    if ($1==chr && $2<=end) {
        if ($3>end) end=$3
    } else {
        print chr, start, end, "segdup", 100, ".", start, end, "'"$COLOR_SEGDUP"'"
        chr=$1; start=$2; end=$3
    }
}
END {print chr, start, end, "segdup", 100, ".", start, end, "'"$COLOR_SEGDUP"'"}
' > segdup.tmp

# 3. Combine all, sort by position
cat gap_regions.tmp segdup.tmp | sort -k1,1V -k2,2n > "${OUTPUT}.tmp"

# 4. Add track header
echo 'track name="hg38_exclude_v1.1" description="hg38 centromere/telomere/heterochromatin/short_arm/segdup regions for VCF filtering" itemRgb="On" visibility="1"' | \
    cat - "${OUTPUT}.tmp" > "$OUTPUT"

# Cleanup
rm -f gap_regions.tmp segdup.tmp "${OUTPUT}.tmp"

# Summary
echo "Created: $OUTPUT"
echo "Regions: $(grep -c -v '^track' "$OUTPUT")"
awk 'NR>1 {sum+=$3-$2} END {printf "Size: %d bp (%.1f Mb, %.2f%%)\n", sum, sum/1000000, sum/3088269832*100}' "$OUTPUT"
echo ""
echo "By category:"
awk 'NR>1 {cat[$4]+=$3-$2; cnt[$4]++} END {for(c in cat) printf "  %s: %d regions, %.1f Mb\n", c, cnt[c], cat[c]/1000000}' "$OUTPUT"

echo "Exit status: $?"
