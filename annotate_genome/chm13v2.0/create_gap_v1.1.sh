#!/bin/bash
# Create chm13v2.0_gap BED files
# v1.1: censat + telomere
# v1.1.1: censat + telomere + segdup (>10kb, >99% identity)
# Output: BED9 format with track header, categories distinguished by name and color
# Usage: ./create_gap_v1.1.sh

set -euxo pipefail

SCRIPT_DIR=$(dirname "$(realpath "$0")")
cd "$SCRIPT_DIR"

# Source files
CENSAT="chm13v2.0_censat_v2.0.bed"
TELOMERE="chm13v2.0_telomere.bed"
SEGDUP_RAW="../../data/genome_stratifications/ucsc_segdup/chm13/sedefSegDups.bed"
OUTPUT_V11="chm13v2.0_gap_v1.1.bed"
OUTPUT_V111="chm13v2.0_gap_v1.1.1.bed"

# Colors (RGB)
COLOR_CENSAT="224,0,0"      # Red
COLOR_TELOMERE="0,0,224"    # Blue
COLOR_SEGDUP="0,128,0"      # Green

# 1. Censat: Extract core satellite regions and merge if within 50kb
#    - Include: hor, hsat, dhor, bsat, gsat, censat, rDNA (exclude: mon, ct)
#    - Merge regions within 50kb gap, label as "censat"
# Output as BED9
awk 'NR>1 && $4 ~ /^(hor|hsat|dhor|bsat|gsat|censat|rDNA)/' "$CENSAT" | \
    cut -f1-3 | sort -k1,1V -k2,2n | awk -v color="$COLOR_CENSAT" '
BEGIN {OFS="\t"; gap=50000}
NR==1 {chr=$1; start=$2; end=$3; next}
{
    if ($1==chr && $2 <= end + gap) {
        if ($3 > end) end = $3
    } else {
        print chr, start, end, "censat", 100, ".", start, end, color
        chr=$1; start=$2; end=$3
    }
}
END {print chr, start, end, "censat", 100, ".", start, end, color}
' > censat.tmp

# 2. Telomere: output as BED9 (no header in source file)
awk -F'\t' '{
    OFS="\t"
    print $1, $2, $3, "telomere", 100, ".", $2, $3, "'"$COLOR_TELOMERE"'"
}' "$TELOMERE" > telomere.tmp

# 3. Segdup: filter by length>=10kb, identity>=99%, merge, output as BED9
awk -F'\t' '
BEGIN {OFS="\t"}
{
    chrom = $1; start = $2; end = $3
    len = end - start
    fracMatch = $24
    if (len >= 10000 && fracMatch >= 0.99) {
        print chrom, start, end
    }
}' "$SEGDUP_RAW" | sort -k1,1V -k2,2n | awk -v color="$COLOR_SEGDUP" '
BEGIN {OFS="\t"}
NR==1 {chr=$1; start=$2; end=$3; next}
{
    if ($1==chr && $2<=end) {
        if ($3>end) end=$3
    } else {
        print chr, start, end, "segdup", 100, ".", start, end, color
        chr=$1; start=$2; end=$3
    }
}
END {print chr, start, end, "segdup", 100, ".", start, end, color}
' > segdup.tmp

# 4. Create v1.1 (censat + telomere only)
cat censat.tmp telomere.tmp | sort -k1,1V -k2,2n > "${OUTPUT_V11}.tmp"
echo 'track name="chm13_gap_v1.1" description="CHM13 censat/telomere regions for VCF filtering" itemRgb="On" visibility="1"' | \
    cat - "${OUTPUT_V11}.tmp" > "$OUTPUT_V11"

# 5. Create v1.1.1 (censat + telomere + segdup)
cat censat.tmp telomere.tmp segdup.tmp | sort -k1,1V -k2,2n > "${OUTPUT_V111}.tmp"
echo 'track name="chm13_gap_v1.1.1" description="CHM13 censat/telomere/segdup regions for VCF filtering" itemRgb="On" visibility="1"' | \
    cat - "${OUTPUT_V111}.tmp" > "$OUTPUT_V111"

# Cleanup
rm -f censat.tmp telomere.tmp segdup.tmp "${OUTPUT_V11}.tmp" "${OUTPUT_V111}.tmp"

# Summary
for OUTPUT in "$OUTPUT_V11" "$OUTPUT_V111"; do
    echo ""
    echo "=== $OUTPUT ==="
    echo "Regions: $(grep -c -v '^track' "$OUTPUT")"
    awk 'NR>1 {sum+=$3-$2} END {printf "Size: %d bp (%.1f Mb, %.2f%%)\n", sum, sum/1000000, sum/3100000000*100}' "$OUTPUT"
    echo "By category:"
    awk 'NR>1 {cat[$4]+=$3-$2; cnt[$4]++} END {for(c in cat) printf "  %s: %d regions, %.1f Mb\n", c, cnt[c], cat[c]/1000000}' "$OUTPUT"
done

echo ""
echo "Exit status: $?"
