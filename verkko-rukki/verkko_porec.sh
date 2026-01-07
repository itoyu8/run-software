#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J verkko_porec
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=30G
#SBATCH -c 16

# Usage: sbatch verkko_porec.sh --hifi <hifi.fastq.gz> --nano <ont.fastq.gz> --porec <porec.fastq.gz> -d <output_dir>
# !!CAUTION!! use -d ~/absolute_path/verkko_output

set -euxo pipefail

# Parse arguments
HIFI=""
NANO=""
POREC=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --hifi)
            HIFI="$2"
            shift 2
            ;;
        --nano)
            NANO="$2"
            shift 2
            ;;
        --porec)
            POREC="$2"
            shift 2
            ;;
        -d|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown option $1"
            echo "Usage: $0 --hifi <hifi.fastq.gz> --nano <ont.fastq.gz> --porec <porec.fastq.gz> -d <output_dir>"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$HIFI" ] || [ -z "$NANO" ] || [ -z "$POREC" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 --hifi <hifi.fastq.gz> --nano <ont.fastq.gz> --porec <porec.fastq.gz> -d <output_dir>"
    exit 1
fi

# Container path
CONTAINER="/home/itoyu8/singularity/verkko-v2.2.1.sif"

# Run Verkko
time singularity exec \
    --bind /home/itoyu8/:/home/itoyu8/ \
    --bind /lustre1:/lustre1/ \
    "$CONTAINER" verkko \
    -d "$OUTPUT_DIR" \
    --hifi "$HIFI" \
    --nano "$NANO" \
    --porec "$POREC" \
    --screen-human-contaminants

echo "Exit status: $?"
