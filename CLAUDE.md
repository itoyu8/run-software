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
- Place output files in same directory as input files
- Use `$(dirname "$INPUT_FILE")` to get input directory
- Use `$(basename "$INPUT_FILE" .ext)` to get base name without extension

## Memory and CPU Settings
Standard SBATCH settings:
- Memory: `--mem-per-cpu=4G`
- CPUs: `-c 32`
- Queue: `-p rjobs,mjobs`