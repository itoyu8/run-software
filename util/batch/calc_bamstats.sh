#!/bin/bash
#SBATCH -J bamstats
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 16

# Usage: ./calc_bamstats.sh <reference_fasta> <base_directory>
# Default reference: /home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta
# Default base directory: /home/itoyu8/project/pelt_expansion/rawdata/castle/quilt

REFERENCE=${1:-/home/itoyu8/database/reference/hg38/v0/Homo_sapiens_assembly38.fasta}
BASE_DIR=${2:-/home/itoyu8/project/pelt_expansion/rawdata/castle/quilt}

echo "Processing BAM files with reference: $REFERENCE"
echo "Base directory: $BASE_DIR"

# Process all BAM files (hardcoded paths based on directory structure)

# HiFi samples
echo "Processing HiFi samples..."
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/hifi/H1437/BL1437.hifi.sorted.bam" > "$BASE_DIR/hifi/H1437/BL1437.hifi.sorted.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/hifi/H2009/BL2009.sorted.bam" > "$BASE_DIR/hifi/H2009/BL2009.sorted.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/hifi/HCC1937/HCC1937BL.sorted.bam" > "$BASE_DIR/hifi/HCC1937/HCC1937BL.sorted.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/hifi/HCC1954/HCC1954BL.sorted.bam" > "$BASE_DIR/hifi/HCC1954/HCC1954BL.sorted.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/hifi/Hs578T/normal.sorted.bam" > "$BASE_DIR/hifi/Hs578T/normal.sorted.stats.txt"

# Illumina samples
echo "Processing Illumina samples..."
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/illumina/H1437/BL1437.markdup.bam" > "$BASE_DIR/illumina/H1437/BL1437.markdup.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/illumina/H2009/BL2009.markdup.bam" > "$BASE_DIR/illumina/H2009/BL2009.markdup.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/illumina/HCC1937/HCC1937BL.markdup.bam" > "$BASE_DIR/illumina/HCC1937/HCC1937BL.markdup.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/illumina/HCC1954/HCC1954BL.markdup.bam" > "$BASE_DIR/illumina/HCC1954/HCC1954BL.markdup.stats.txt"

# ONT samples
echo "Processing ONT samples..."
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/ont/H1437/BL1437.ont.sorted.bam" > "$BASE_DIR/ont/H1437/BL1437.ont.sorted.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/ont/H2009/BL2009.ont.sorted.bam" > "$BASE_DIR/ont/H2009/BL2009.ont.sorted.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/ont/HCC1395/HCC1395BL.sorted.ont.bam" > "$BASE_DIR/ont/HCC1395/HCC1395BL.sorted.ont.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/ont/HCC1937/HCC1937BL.sorted.ont.bam" > "$BASE_DIR/ont/HCC1937/HCC1937BL.sorted.ont.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/ont/HCC1954/HCC1954BL.sorted.ont.bam" > "$BASE_DIR/ont/HCC1954/HCC1954BL.sorted.ont.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/ont/Hs578T/Hs578Bst.ont.sorted.bam" > "$BASE_DIR/ont/Hs578T/Hs578Bst.ont.sorted.stats.txt"

# UL samples
echo "Processing UL samples..."
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/ul/H1437/BL1437.ul.sorted.bam" > "$BASE_DIR/ul/H1437/BL1437.ul.sorted.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/ul/H2009/BL2009.ul.sorted.bam" > "$BASE_DIR/ul/H2009/BL2009.ul.sorted.stats.txt"
/home/itoyu8/bin/samtools/samtools-1.19/samtools stats -@ 16 --reference "$REFERENCE" "$BASE_DIR/ul/HCC1937/HCC1937BL.ul.sorted.bam" > "$BASE_DIR/ul/HCC1937/HCC1937BL.ul.sorted.stats.txt"

echo "All BAM stats calculations completed."