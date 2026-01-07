#!/bin/bash
# Calculate remaining bases after excluding gap regions
#
# === Calculation Logic ===
# CHM13: remaining = total - gap_v1.1(censat + telomere)
#        - CHM13 is fully sequenced (no Ns)
#        - gap_v1.1 contains censat (50kb merged) + telomere
#
# hg38:  remaining = total - union(N_regions, gap_v1.0)
#        - N_regions: from UCSC gap.txt (contig, scaffold, heterochromatin, short_arm, telomere)
#        - gap_v1.0: centromere, heterochromatin, short_arm, telomere
#        - Note: centromere is NOT N-masked in hg38 (it's model sequence)
#        - Note: Some contig/scaffold gaps fall within centromere (269,500 bp)
#        - overlap = heterochromatin + short_arm + telomere + (contig/scaffold in centromere)
#                  = 137,377,000 + 269,500 = 137,646,500 bp
#
# === Output ===
# Per-chromosome breakdown for chr1-22, chrX, chrY (alt contigs excluded)
#
# Usage: ./calc_remaining.sh

set -uo pipefail

SCRIPT_DIR=$(dirname "$(realpath "$0")")

# File paths
CHM13_SIZES="${SCRIPT_DIR}/chm13v2.0/chm13v2.0.chrom.sizes"
CHM13_GAP="${SCRIPT_DIR}/chm13v2.0/chm13v2.0_gap_v1.1.bed"
HG38_SIZES="${SCRIPT_DIR}/hg38/hg38.chrom.sizes"
HG38_N_REGIONS="${SCRIPT_DIR}/hg38/hg38_gap_ucsc.txt"
HG38_GAP="${SCRIPT_DIR}/hg38/hg38_gap_v1.0.bed"

# Temp files
TMP_DIR=$(mktemp -d)
trap "rm -rf $TMP_DIR" EXIT

###############################################################################
# CHM13: remaining = total - gap(censat + telomere)
###############################################################################
echo "=== CHM13 v1.1 ==="
echo "計算式: remaining = total - gap(censat + telomere)"
echo ""
echo -e "chr\ttotal\tgap\tremaining\t%"

awk '
BEGIN {OFS="\t"}
NR==FNR {
    size[$1] = $2
    next
}
$4=="censat" || $4=="telomere" {
    gap[$1] += $3 - $2
}
END {
    total_size = 0
    total_gap = 0
    for (chr in size) {
        if (chr ~ /^chr([0-9]+|X|Y)$/) {
            remaining = size[chr] - gap[chr]
            pct = remaining / size[chr] * 100
            printf "%s\t%d\t%d\t%d\t%.2f\n", chr, size[chr], gap[chr], remaining, pct
            total_size += size[chr]
            total_gap += gap[chr]
        }
    }
    remaining = total_size - total_gap
    pct = remaining / total_size * 100
    printf "TOTAL\t%d\t%d\t%d\t%.2f\n", total_size, total_gap, remaining, pct
}
' "$CHM13_SIZES" "$CHM13_GAP" | sort -k1,1V > "${TMP_DIR}/chm13_result.txt"
grep -v "^TOTAL" "${TMP_DIR}/chm13_result.txt"
grep "^TOTAL" "${TMP_DIR}/chm13_result.txt"

###############################################################################
# hg38: remaining = total - union(N_regions, gap_v1.0)
###############################################################################
echo ""
echo "=== hg38 v1.0 ==="
echo "計算式: remaining = total - union(N_regions, gap_v1.0)"
echo ""

# Step 1: Combine N regions and gap regions, then merge overlaps
# This correctly handles:
# - Overlapping regions (e.g., heterochromatin in both N and gap)
# - contig/scaffold gaps within centromere regions
{
    # N regions from UCSC gap.txt (columns: bin, chrom, start, end, ...)
    awk '$2 ~ /^chr([0-9]+|X|Y)$/ {print $2"\t"$3"\t"$4}' "$HG38_N_REGIONS"
    # gap_v1.0 regions (skip header)
    awk 'NR>1 && $1 ~ /^chr([0-9]+|X|Y)$/ {print $1"\t"$2"\t"$3}' "$HG38_GAP"
} | sort -k1,1V -k2,2n | awk '
BEGIN {OFS="\t"}
{
    if (NR == 1) {
        chr = $1; start = $2; end = $3
    } else if ($1 == chr && $2 <= end) {
        # Overlapping or adjacent - extend
        if ($3 > end) end = $3
    } else {
        # New region - output previous
        print chr, start, end
        chr = $1; start = $2; end = $3
    }
}
END {
    print chr, start, end
}
' > "${TMP_DIR}/hg38_merged.bed"

# Step 2: Calculate union size per chromosome
awk '
BEGIN {OFS="\t"}
{
    merged[$1] += $3 - $2
}
END {
    for (chr in merged) {
        print chr, merged[chr]
    }
}
' "${TMP_DIR}/hg38_merged.bed" > "${TMP_DIR}/hg38_union.txt"

# Step 3: Calculate N regions per chromosome
awk '
BEGIN {OFS="\t"}
$2 ~ /^chr([0-9]+|X|Y)$/ {
    n[$2] += $4 - $3
}
END {
    for (chr in n) {
        print chr, n[chr]
    }
}
' "$HG38_N_REGIONS" > "${TMP_DIR}/hg38_n.txt"

# Step 4: Calculate gap_v1.0 per chromosome
awk '
BEGIN {OFS="\t"}
NR > 1 && $1 ~ /^chr([0-9]+|X|Y)$/ {
    gap[$1] += $3 - $2
}
END {
    for (chr in gap) {
        print chr, gap[chr]
    }
}
' "$HG38_GAP" > "${TMP_DIR}/hg38_gap.txt"

# Step 5: Output results
echo -e "chr\ttotal\tN_regions\tgap_v1.0\toverlap\tunion\tremaining\t%"

awk '
BEGIN {OFS="\t"}
FILENAME ~ /chrom.sizes/ {
    size[$1] = $2
    next
}
FILENAME ~ /hg38_n.txt/ {
    n[$1] = $2
    next
}
FILENAME ~ /hg38_gap.txt/ {
    gap[$1] = $2
    next
}
FILENAME ~ /hg38_union.txt/ {
    merged[$1] = $2
    next
}
END {
    total_size = 0
    total_n = 0
    total_gap = 0
    total_merged = 0

    for (chr in size) {
        if (chr ~ /^chr([0-9]+|X|Y)$/) {
            n_val = (chr in n) ? n[chr] : 0
            gap_val = (chr in gap) ? gap[chr] : 0
            merged_val = (chr in merged) ? merged[chr] : 0
            overlap = n_val + gap_val - merged_val
            remaining = size[chr] - merged_val
            pct = remaining / size[chr] * 100
            printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", chr, size[chr], n_val, gap_val, overlap, merged_val, remaining, pct

            total_size += size[chr]
            total_n += n_val
            total_gap += gap_val
            total_merged += merged_val
        }
    }

    total_overlap = total_n + total_gap - total_merged
    total_remaining = total_size - total_merged
    total_pct = total_remaining / total_size * 100
    printf "TOTAL\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", total_size, total_n, total_gap, total_overlap, total_merged, total_remaining, total_pct
}
' "$HG38_SIZES" "${TMP_DIR}/hg38_n.txt" "${TMP_DIR}/hg38_gap.txt" "${TMP_DIR}/hg38_union.txt" | sort -k1,1V > "${TMP_DIR}/hg38_result.txt"
grep -v "^TOTAL" "${TMP_DIR}/hg38_result.txt"
grep "^TOTAL" "${TMP_DIR}/hg38_result.txt"

###############################################################################
# Summary
###############################################################################
echo ""
echo "=== Summary ==="
echo "CHM13: remaining = total - gap(censat + telomere)"
echo "hg38:  remaining = total - union(N_regions, gap_v1.0)"
echo ""
echo "hg38 overlap breakdown:"
echo "  - Same categories (heterochromatin + short_arm + telomere): 137,377,000 bp"
echo "  - contig/scaffold within centromere: 269,500 bp"
echo "  - Total overlap: 137,646,500 bp"
