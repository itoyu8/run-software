# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Repository Overview

This is a bioinformatics pipeline repository containing SLURM job scripts for genomic analysis workflows. The repository is organized into tool-specific directories, each containing scripts for different stages of genomic data processing.

## Directory Structure

- **bwa-samtools/**: Short-read and long-read alignment pipelines using BWA/Minimap2 and Samtools processing
- **glimpse1/**, **glimpse2/**: Genotype imputation using GLIMPSE (versions 1 and 2)
- **quilt/**: Genotype imputation using QUILT2
- **rasusa/**: FASTQ downsampling using Rasusa
- **seqtk/**: FASTQ downsampling using seqtk
- **dv-whatshap/**: DeepVariant and WhatsHap phasing pipelines
- **svcaller/**: Structural variant calling (Severus, nanomonsv)
- **gatk/**: GATK utility scripts
- **temp/**: Work-in-progress or experimental scripts

## Key Workflow Patterns

### Imputation Workflows
Both GLIMPSE2 and QUILT follow a three-stage pattern:
1. **Reference preparation**: Convert 1000 Genomes VCF to tool-specific format
2. **Chunking**: Divide genome into manageable regions using genetic maps
3. **Imputation**: Process BAM files per chromosome/chunk, then ligate/merge results

**GLIMPSE2 execution order:**
```bash
sbatch prepare_refpanel.sh  # One-time setup
sbatch make_chunks.sh       # One-time setup
sbatch split_reference.sh   # One-time setup
sbatch run_glimpse2.sh /path/to/sample.bam [output_name]
```

**QUILT execution order:**
```bash
sbatch create_chunks.sh       # One-time setup
sbatch prepare_reference.sh   # One-time setup
sbatch run_quilt.sh /path/to/sample.bam [output_name]
```

### Alignment Workflows
- **Short-reads**: `bwa_samsort.sh` performs BWA alignment → Samtools sort → GATK MarkDuplicates
- **Long-reads**: `minimap2_samsort.sh` supports ONT and HiFi with `--type` parameter
- Both scripts output BAM, BAI, and optionally metrics files in the same directory as input

### Container Usage Strategy
- **Prefer direct binaries** when available (BWA, Samtools, Minimap2, BCFtools)
- **Use Singularity containers** for tools without local binaries (GLIMPSE, QUILT, GATK)
- Container paths are hard-coded in scripts (e.g., `/home/itoyu8/singularity/glimpse_v2.0.0-27-g0919952_20221207.sif`)

## Common Script Parameters

### Reference Genome Selection
Most scripts accept `--reference hg38|chm13` to switch between GRCh38 and CHM13 reference genomes.

### Output Conventions
- Output files are placed in the **same directory as input files** by default
- Many scripts accept an optional output name/folder parameter (e.g., `run_glimpse2.sh`, `run_quilt.sh`)
- Output names can include subdirectories (e.g., `results/analysis_name`)

### Downsampling Tools
- **rasusa**: Coverage-based downsampling (specify target coverage like `--coverage 30`)
- **seqtk**: Fraction-based downsampling (specify sampling rate like `-r 0.1`)
- **samtools**: BAM downsampling via `sam_downsample.sh`

## Singularity Container Locations

Key containers used in scripts:
- GLIMPSE2: `/home/itoyu8/singularity/glimpse_v2.0.0-27-g0919952_20221207.sif`
- QUILT: `/home/itoyu8/singularity/quilt_v0.1.0.sif`
- GATK: `/home/itoyu8/singularity/compat_parabricks-0.2.2.sif`

## Testing and Validation

When modifying scripts:
1. Verify SLURM directives are intact (especially memory and CPU requirements)
2. Check that input file path handling preserves directory structure
3. Ensure output files are created in the correct location (typically same directory as input)
4. Test with both hg38 and chm13 reference options if applicable
5. Verify log files will be created in `./log/` directory (must exist before job submission)

# Coding Standards for Bioinformatics Scripts

## Shebang Format
Always use SBATCH directives in the shebang:
```bash
#!/bin/bash
#SBATCH -p rjobs,mjobs
#SBATCH -J script_name
#SBATCH -o ./log/%x.o%j
#SBATCH -e ./log/%x.e%j
#SBATCH --mem-per-cpu=4G
#SBATCH -c 32
```

## Usage Documentation
- Use `# Usage:` comment at line 8 instead of extensive error handling
- Keep usage examples concise and clear
- Example: `# Usage: ./script.sh [--reference hg38|chm13] <input_file>`

## Reference Genome Paths
Standard reference genome paths and CHM13 branching:
```bash
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.fa"
else
    REFERENCE_GENOME_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"
fi
```

For minimap2 MMI files:
```bash
if [ "$REFERENCE_TYPE" = "chm13" ]; then
    REFERENCE_MMI_PATH="/home/itoyu8/database/reference/chm13/v2.0/chm13v2.0_maskedY_rCRS.mmi"
else
    REFERENCE_MMI_PATH="/home/itoyu8/database/reference/hg38/GRCh38.d1.vd1/GRCh38.d1.vd1.mmi"
fi
```

## Thread Definition
Always use SLURM variable with fallback:
```bash
THREADS=${SLURM_CPUS_PER_TASK:-32}
```

## Command Paths

### Direct Binary Paths (Preferred)
- BWA: `/home/itoyu8/bin/bwa/bwa-0.7.19/bwa`
- Samtools: `/home/itoyu8/bin/samtools/samtools-1.19/samtools`
- Minimap2: `/home/itoyu8/bin/minimap2/minimap2-2.28/minimap2`

### Container Usage
Only use containers when direct binaries are not available:
- GATK container: `/home/itoyu8/singularity/compat_parabricks-0.2.2.sif`
- Container execution: `singularity exec ${CONTAINER_PATH} gatk`

## Argument Parsing
Use consistent argument parsing pattern:
```bash
INPUT_FILES=()
REFERENCE_TYPE="hg38"

while [[ $# -gt 0 ]]; do
    case $1 in
        --reference)
            if [ "$2" = "chm13" ]; then
                REFERENCE_TYPE="chm13"
            elif [ "$2" = "hg38" ]; then
                REFERENCE_TYPE="hg38"
            else
                echo "Error: --reference must be 'hg38' or 'chm13'"
                exit 1
            fi
            shift 2
            ;;
        *)
            INPUT_FILES+=("$1")
            shift
            ;;
    esac
done
```

## File Output

### CRITICAL: Output Directory Handling
**When a script accepts an output path parameter (e.g., `-o`, `-d`), treat it as a DIRECTORY, not a file prefix.**

Common mistake pattern to AVOID:
```bash
# WRONG - This places files in parent directory
OUTPUT_DIR=$(dirname "${OUTPUT_PREFIX}")
mkdir -p "${OUTPUT_DIR}"
command -o "${OUTPUT_PREFIX}" ...  # Files end up scattered in parent directory
```

Correct pattern:
```bash
# CORRECT - Treat output parameter as directory
OUTPUT_DIR="${OUTPUT_PREFIX}"  # The parameter IS the directory
mkdir -p "${OUTPUT_DIR}"
OUTPUT_BASE=$(basename "${OUTPUT_PREFIX}")
command -o "${OUTPUT_DIR}/${OUTPUT_BASE}" ...  # All files inside the directory
```

**Examples of tools following this pattern:**
- **hifiasm**: `-o ~/path/output_name` → All files in `output_name/` directory
- **dipcall**: `-o ~/path/output_name` → All files in `output_name/` directory
- **verkko**: `-d ~/path/output_dir` → All files in `output_dir/` directory

**Why this matters:**
- Users expect `-o sample1` to create `sample1/` containing all outputs
- The alternative (files scattered in parent directory) causes organization issues
- Consistent with assembly tool conventions (hifiasm, verkko, etc.)

### Legacy Pattern (for alignment scripts only)
- Place output files in same directory as input files
- Use `$(dirname "$INPUT_FILE")` to get input directory
- Use `$(basename "$INPUT_FILE" .ext)` to get base name without extension

## Memory and CPU Settings
Standard SBATCH settings:
- Memory: `--mem-per-cpu=4G`
- CPUs: `-c 32`
- Queue: `-p rjobs,mjobs`

## Docker Image Management

### Building and Pushing to Docker Hub
IMPORTANT: Always use a two-step tagging process. Direct tagging may not be recognized properly.

```bash
# Step 1: Build with platform specification and temporary tag
docker build --platform linux/amd64 -t [image_name]-amd64:[version] -f ./Dockerfile .

# Step 2: Re-tag for Docker Hub (this step is REQUIRED, don't skip)
# Direct push without re-tagging may fail to be recognized
docker tag [image_name]-amd64:[version] itoyu8/[image_name]:[version]

# Step 3: Push to Docker Hub
docker push itoyu8/[image_name]:[version]
```

Example workflow:
```bash
# Build rasusa image
docker build --platform linux/amd64 -t rasusa-amd64:0.1.0 -f ./Dockerfile .

# Re-tag (REQUIRED - don't skip this intermediate step)
docker tag rasusa-amd64:0.1.0 itoyu8/rasusa:0.1.0

# Push to Docker Hub
docker push itoyu8/rasusa:0.1.0
```

### Why Two-Step Tagging?
- Past experience shows that skipping the intermediate tagging step causes recognition issues
- The re-tagging step ensures Docker properly registers the image for push
- Always go through the temporary tag → Docker Hub tag workflow

## File Transfer with SCP

When requested by the user under their supervision, use `scp` to transfer files to remote servers:

```bash
scp <local_file> itoyu8@10.20.27.13:/home/itoyu8/project/pelt_expansion/rscarpia_validate/
```

**Important Notes:**
- Only use scp when the user explicitly requests it and is actively monitoring
- The user has configured SSH public key authentication for passwordless access
- Common transfer targets:
  - Development server: `itoyu8@10.20.27.13:/home/itoyu8/project/pelt_expansion/rscarpia_validate/`