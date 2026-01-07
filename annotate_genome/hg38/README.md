# hg38 Exclude Regions for VCF Filtering

## Overview

This directory contains BED files defining genomic regions that should be excluded from variant analysis due to low mappability or unreliable variant calling. These regions are specifically designed for use with the **GRCh38.d1.vd1** reference genome (GDC reference).

## Source Data

All annotation files were downloaded from UCSC Genome Browser database (hg38):
- **URL**: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/
- **Download date**: 2024-12-11

### Original Files

| File | Description | Source |
|------|-------------|--------|
| `centromeres.txt` | Centromere model sequence positions | UCSC centromeres track |
| `gap.txt` | Assembly gap annotations | UCSC gap track |

## Processing Methods

### Centromere Regions

The hg38 assembly contains **modeled centromere sequences** (alpha-satellite representations) rather than actual sequenced centromeres. The `centromeres.txt` file contains multiple contigs per chromosome representing these model sequences.

**Processing steps:**
1. Extracted chromosome, start, and end positions from `centromeres.txt`
2. For each chromosome, calculated the **minimum start** and **maximum end** positions to create a single continuous region that encompasses all centromere model contigs
3. This approach ensures complete coverage of the centromeric region, including small gaps between model sequence contigs

```bash
awk 'BEGIN{OFS="\t"}
{
  if (!($1 in min) || $2 < min[$1]) min[$1] = $2
  if (!($1 in max) || $3 > max[$1]) max[$1] = $3
}
END {
  for (chr in min) print chr, min[chr], max[chr], "centromere"
}' centromeres.txt | sort -k1,1V > hg38_centromeres_minmax.bed
```

### Gap Regions

From `gap.txt`, the following gap types were extracted for exclusion:

| Gap Type | Description | Reason for Exclusion |
|----------|-------------|---------------------|
| `telomere` | Chromosome ends | Repetitive sequences, unreliable mapping |
| `heterochromatin` | Constitutive heterochromatin (e.g., 1q12, 9q12, 16q11.2) | Highly repetitive, poor mappability |
| `short_arm` | Acrocentric chromosome short arms (chr13, 14, 15, 21, 22 p-arms) | rDNA repeats, no unique sequence |

```bash
awk '$8=="telomere" || $8=="heterochromatin" || $8=="short_arm" {
  print $2, $3, $4, $8
}' gap.txt > hg38_gap_exclude.bed
```

### Final Exclude Regions

The final file combines all exclude regions:

```bash
cat hg38_centromeres_minmax.bed hg38_gap_exclude.bed | sort -k1,1V -k2,2n > hg38_exclude_regions_final.bed
```

## Output Files

| File | Description | Format |
|------|-------------|--------|
| `hg38_censat_gap_v1.0.bed` | **Combined exclude regions for VCF filtering** | BED9 with track header |

### Intermediate Files (in `intermediate/` directory)

| File | Description | Format |
|------|-------------|--------|
| `centromeres.txt` | Original UCSC centromeres data | UCSC table |
| `gap.txt` | Original UCSC gap data | UCSC table |
| `hg38_centromeres.bed` | Raw centromere contigs from UCSC | BED9 |
| `hg38_centromeres_minmax.bed` | Merged centromere regions (min-max per chromosome) | BED9 |
| `hg38_gap_exclude.bed` | Telomere, heterochromatin, and short_arm regions | BED9 |
| `hg38_exclude_regions.bed` | Combined regions without header | BED9 |

## Region Summary

| Region Type | Count | Description |
|-------------|-------|-------------|
| centromere | 24 | One per chromosome (chr1-22, X, Y) |
| telomere | 48 | Two per chromosome (p-arm and q-arm ends) |
| heterochromatin | 11 | Large heterochromatic blocks |
| short_arm | 5 | Acrocentric chromosome short arms |
| **Total** | **88** | |

## Usage

For VCF filtering with bcftools:
```bash
bcftools view -T ^hg38_censat_gap_v1.0.bed input.vcf.gz -Oz -o filtered.vcf.gz
```

For bedtools intersection:
```bash
bedtools intersect -v -a input.vcf -b hg38_censat_gap_v1.0.bed > filtered.vcf
```

## Important Notes

1. **Coordinate system**: All coordinates are 0-based, half-open (BED format), compatible with hg38/GRCh38.
2. **Chromosome naming**: Uses UCSC-style naming with `chr` prefix (chr1, chr2, ..., chrX, chrY, chrM).
3. **Reference compatibility**: Validated for use with `GRCh38.d1.vd1.fa` (GDC reference). Primary chromosome coordinates are identical to standard UCSC hg38.

## References

- UCSC Genome Browser: https://genome.ucsc.edu/
- GRCh38 Assembly: GCA_000001405.15
- GDC Reference: https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
