# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Repository Overview

Bioinformatics pipeline repository containing SLURM job scripts for genomic analysis workflows on HPC. Each tool-specific directory contains scripts for different stages of genomic data processing.

Coding standards (shebang, shell options, exit status, echo rules, time measurement, reference paths, thread definition, command paths, argument parsing, output conventions, SBATCH settings, Docker, SCP) are defined in the global CLAUDE.md and apply to all scripts here.

## Directory Organization

- **Alignment**: `bwa-samtools/` (BWA, Minimap2, samtools stats/downsample/bedcov), `bam_refiner/`
- **Germline variant calling**: `dv-whatshap/` (DeepVariant + WhatsHap), `dvr9-whatshap/` (DV R9)
- **Somatic variant calling**: `clairs/` (ClairS), `deepsomatic/` (DeepSomatic), `svcaller/` (Severus, nanomonsv), `wakhan/` (somatic CNA profiling)
- **Imputation**: `glimpse1/`, `glimpse2/`, `quilt/hg38/`, `quilt/chm13/`, `beagle/`
- **Assembly**: `hifiasm/`, `verkko-rukki/`, `dipcall/`, `switch_error/`
- **Downsampling**: `rasusa/` (coverage-based), `seqtk/` (fraction-based)
- **Utilities**: `gatk/`, `mocha/`, `annotate_genome/`, `process_snp/`, `util/`

## Key Workflow Patterns

### Alignment

- **Short-reads**: `bwa-samtools/bwa_samsort.sh` (BWA align -> samtools sort -> GATK MarkDuplicates)
- **Long-reads**: `bwa-samtools/minimap2_samsort.sh` supports ONT and HiFi via `--type` parameter

### Germline Variant Calling + Phasing (DeepVariant + WhatsHap)

```bash
# Step 1: Variant calling and phasing
sbatch dv-whatshap/dv_whphase.sh --type ont -d /output/dir -o sample_name /path/to/sample.bam
# -> sample_name.dv.vcf.gz, sample_name.phased.vcf.gz

# Step 2: Haplotagging and BAM splitting
sbatch dv-whatshap/whtag_split.sh -d /output/dir sample_name.phased.vcf.gz sample_name.bam
# -> sample_name.hptag.bam, sample_name.h1.bam, sample_name.h2.bam
```

### DeepVariant Rescue (CHM13 only)

CHM13 telomere regions can cause DeepVariant "invalid allele index" crashes, losing data for the affected chromosome and subsequent chromosomes in the same shard. The rescue workflow (`dv-whatshap/dv_rescue.sh`) re-runs DeepVariant on affected regions and merges results. See `dv-whatshap/DV_RESCUE_GUIDE.md` for the full diagnostic and recovery procedure.

```bash
bash dv-whatshap/dv_rescue.sh \
    --regions "chr17:30,chr18,chr21,chr22" \
    --original-vcf <path>/normal.dv.vcf.gz \
    --type ont --reference chm13 --strict-filter \
    -d <path>/rescue <path>/normal.bam
```

### Somatic Variant Calling (tumor/normal pairs)

```bash
# ClairS (with optional pre-phased VCF)
bash clairs/clairs.sh -d /output/dir tumor.bam normal.bam
bash clairs/clairs.sh -d /output/dir --normal-vcf germline.vcf.gz --haplotagged tumor_hptag.bam normal.bam

# DeepSomatic
bash deepsomatic/deepsomatic.sh -d /output/dir --platform ont tumor.bam normal.bam

# Severus (structural variants)
./svcaller/severus/severus.sh --tumor tumor.bam --normal normal.bam --phased-vcf phased.vcf --out-dir /output/dir

# Wakhan (somatic CNA profiling, with Severus breakpoints)
sbatch wakhan/wakhan.sh -d /out -o sample --phased-vcf normal.phased.vcf.gz --breakpoints severus_somatic.vcf tumor.bam

# Wakhan (standalone with change-point detection)
sbatch wakhan/wakhan.sh -d /out -o sample --phased-vcf normal.phased.vcf.gz --cpd tumor.bam
```

### Imputation Workflows

GLIMPSE2 and QUILT both follow a three-stage pattern: reference preparation (one-time) -> chunking (one-time) -> per-sample imputation.

**GLIMPSE2:**
```bash
sbatch glimpse2/prepare_refpanel.sh   # One-time
sbatch glimpse2/make_chunks.sh        # One-time
sbatch glimpse2/split_reference.sh    # One-time
sbatch glimpse2/run_glimpse2.sh /path/to/sample.bam [output_name]
```

**QUILT2 (hg38)** - uses SHAPEIT2 phased 1KGP panel (3202 samples):
```bash
sbatch quilt/hg38/create_chunks.sh       # One-time
sbatch quilt/hg38/prepare_reference.sh   # One-time
sbatch quilt/hg38/quilt_multi.sh [-d output_dir] <input.bam>
```

**QUILT2 (chm13)** - uses SHAPEIT5 phased 1KGP panel (2504 unrelated) with T2T-native genetic maps:
```bash
bash quilt/chm13/download_2504_bcf.sh                   # Mac (HPC has no internet)
sbatch quilt/chm13/create_chunks.sh                      # One-time
sbatch quilt/chm13/split_bcf.sh                          # One-time
for i in {1..22}; do sbatch quilt/chm13/prepare_reference.sh chr$i; done  # One-time, parallel
```

### Assembly (hifiasm)

```bash
# Trio mode with yak DBs
sbatch hifiasm/hifiasm_trio_yak.sh --paternal pat.yak --maternal mat.yak -d /out hifi.fastq.gz

# ONT-only mode
sbatch hifiasm/hifiasm_ontonly.sh -d /out ont.fastq.gz

# GFA -> FASTA post-processing
bash hifiasm/gfaprocess.sh <assembly_prefix>
```

## Container Usage Strategy

- **Prefer direct binaries** when available (BWA, Samtools, Minimap2, BCFtools, Yak)
- **Use Singularity containers** for tools without local binaries:
  - DeepVariant+WhatsHap: `dv-whatshap_0.1.0.sif`
  - ClairS: `clairs_0.1.0.sif`
  - DeepSomatic: `deepsomatic_0.1.0.sif`
  - Severus: `severus_0.1.0.sif`
  - Wakhan: `wakhan_0.4.2.sif`
  - GLIMPSE2: `glimpse_v2.0.0-27-g0919952_20221207.sif`
  - QUILT: `quilt_v0.1.0.sif`
  - GATK: `compat_parabricks-0.2.2.sif`
  - Python3: `python3_0.1.0.sif`
- All containers stored at `/home/itoyu8/singularity/`

## Common Script Parameters

- `--reference hg38|chm13` - switch reference genome (default: hg38)
- `--type ont|hifi` - sequencing platform (for alignment and variant calling)
- `-d <dir>` - output directory
- `-o <name>` - output name/prefix

## Testing and Validation

When modifying scripts:
1. Verify SLURM directives are intact (especially memory and CPU requirements)
2. Check that input file path handling preserves directory structure
3. Ensure output files are created in the correct location
4. Test with both hg38 and chm13 reference options if applicable
5. Verify `./log/` directory exists before job submission (logs go to `./log/%x.o%j`)
