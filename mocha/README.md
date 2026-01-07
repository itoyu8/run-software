# MoChA (MOsaic CHromosomal Alterations) Pipeline

Mosaic chromosomal alteration detection from WGS data using MoChA.

## Prerequisites

### Software
- BCFtools 1.19+
- Singularity
- MoChA container (`mocha_0.1.0.sif`)

### Resource Files

Download from Broad Institute:
```bash
wget https://software.broadinstitute.org/software/mocha/mocha.GRCh38.zip
unzip mocha.GRCh38.zip
```

Required files:
- `segdups.bed.gz` - Segmental duplications with Jukes-Cantor divergence
- `cnps.bed` - Known copy number polymorphisms

Resource location: `/home/itoyu8/database/tools/mocha/`

## Prepare Exclusion BED (One-time Setup)

Extract low-divergence segmental duplication regions (JK < 0.02) to exclude from analysis:

```bash
zcat /home/itoyu8/database/tools/mocha/segdups.bed.gz | \
    awk -F'\t' '$4 < 0.02 {print $1"\t"$2"\t"$3}' | \
    bgzip > /home/itoyu8/database/tools/mocha/segdups_exclude.bed.gz

tabix -p bed /home/itoyu8/database/tools/mocha/segdups_exclude.bed.gz
```

This creates a BED file of regions where segmental duplications have very low sequence divergence (< 2%), which can cause false positive mosaic calls.

## Usage

```bash
sbatch mocha.sh <GATK_VCF> <PHASED_VCF> <OUTPUT_PREFIX>
```

### Input Files
- `GATK_VCF`: VCF from GATK with AD (allelic depth) field
- `PHASED_VCF`: Phased VCF (e.g., from Beagle)
- `OUTPUT_PREFIX`: Output directory path

### Output Files
- `mocha_input.vcf.gz` - Input VCF with AD and GC annotations
- `mocha_calls.tsv` - Detected mosaic chromosomal alterations
- `mocha_stats.tsv` - Per-sample genome-wide statistics

## Pipeline Steps

1. **AD Annotation**: Transfer AD field from GATK VCF to phased VCF
2. **GC Annotation**: Add GC content field using `bcftools +mochatools`
3. **MoChA Analysis**: Detect mosaic chromosomal alterations

## MoChA Options Used

| Option | Description |
|--------|-------------|
| `-g GRCh38` | Use GRCh38 genome rules (centromeres, MHC, KIR regions) |
| `-p cnps.bed` | Known CNP regions (flag overlapping calls) |
| `-T ^segdups_exclude.bed.gz` | Exclude low-divergence segdup regions |
| `-z stats.tsv` | Output per-sample statistics |
| `-c calls.tsv` | Output mosaic calls |

## References

- MoChA: https://github.com/freeseek/mocha
- Loh PR, et al. Insights into clonal haematopoiesis from 8,342 mosaic chromosomal alterations. Nature. 2018.
