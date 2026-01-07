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

set -euxo pipefail
```

## Shell Options (TODO: Apply to all scripts)
All scripts should include `set -euxo pipefail` after SBATCH directives:
- `-e`: Exit on error
- `-u`: Error on undefined variables
- `-x`: Print commands (for debugging which sample failed)
- `-o pipefail`: Catch errors in pipes

Note: These options are NOT inherited by child scripts when using `bash script.sh`. Each script needs its own `set` line.

## Exit Status
Always print the exit status at the end of scripts for error log identification:
```bash
echo "Exit status: $?"
```

This helps identify failed jobs when reviewing logs. With `set -e`, if a command fails the script exits immediately, so reaching this line indicates success (exit status 0).

## Echo Statements
With `set -x` enabled, all commands are printed with expanded variables before execution. Therefore, explicit `echo` statements for debugging or showing variable values are unnecessary and should be avoided. Only use `echo` for:
- Error messages before `exit 1`
- Usage information
- Exit status at the end of scripts

## Time Measurement
Use the `time` command to measure execution time instead of manual `START_TIME`/`END_TIME` calculations:
```bash
# Good - use time command
time "${TOOL}" args > output.log 2>&1

# Bad - manual time calculation (deprecated)
START_TIME=$(date +%s)
"${TOOL}" args
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
```

The `time` command outputs:
- `real`: Wall clock time (actual elapsed time)
- `user`: CPU time in user mode (sum across all threads)
- `sys`: CPU time in kernel mode

For multi-threaded tools, `user` time will be greater than `real` time, which indicates parallel efficiency.

## Usage Documentation
- Use `# Usage:` comment at line 8 instead of extensive error handling
- Use `# Output:` comment at line 9 to document output files
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
- Yak: `/home/itoyu8/bin/yak/yak-0.1/yak`
- k8 (for paftools.js): `/home/itoyu8/bin/minimap2/k8-0.2.5/k8-Linux`

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

### Output Path Convention (TODO: Standardize across all scripts)
Use `-d` for output directory and `-o` for output name/prefix:
```bash
# Usage pattern
script.sh -d /output/dir -o sample_name input.bam
# → /output/dir/sample_name.vcf.gz

# -o has a default, so -d alone works
script.sh -d /output/dir input.bam
# → /output/dir/output.vcf.gz (default name)
```

Implementation pattern:
```bash
OUTPUT_DIR="."
OUTPUT_NAME="output"

while [[ $# -gt 0 ]]; do
    case $1 in
        -d) OUTPUT_DIR="$2"; shift 2 ;;
        -o) OUTPUT_NAME="$2"; shift 2 ;;
        *)  INPUT_FILE="$1"; shift ;;
    esac
done

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(realpath "${OUTPUT_DIR}")
# → Use ${OUTPUT_DIR}/${OUTPUT_NAME}.ext for outputs
```

### Current issues (TODO)
1. Some scripts use realpath, others don't → Standardize to always use realpath
2. Output directory not auto-created → Always `mkdir -p`
3. Mix of output dir vs output filename → Migrate to `-d` / `-o` pattern above

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