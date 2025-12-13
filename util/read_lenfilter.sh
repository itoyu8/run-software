#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J read_lenfilter
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=8G
#SBATCH -c 1
# Usage: sbatch read_lenfilter.sh <INPUT_FASTQ> -m <MIN_LENGTH> -o <OUTPUT_FASTQ>

INPUT_FASTQ=$1
MIN_LENGTH=$3
OUTPUT_FASTQ=$5

mkdir -p log

/home/itoyu8/bin/seqkit/seqkit seq -m "${MIN_LENGTH}" "${INPUT_FASTQ}" > "${OUTPUT_FASTQ}"
