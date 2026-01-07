#!/usr/bin/env python3
"""
Plot contig alignments from PAF file for a specific chromosome.
Usage: python plot_paf.py --reference chm13|hg38 --chr chrN input.paf
Output: input.chrN.png (in same directory as input PAF)
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# Chromosome lengths for references
CHM13_CHROM_LENGTHS = {
    "chr1": 248387328, "chr2": 242696752, "chr3": 201105948, "chr4": 193574945,
    "chr5": 182045439, "chr6": 172126628, "chr7": 160567428, "chr8": 146259331,
    "chr9": 150617247, "chr10": 134758134, "chr11": 135127769, "chr12": 133324548,
    "chr13": 113566686, "chr14": 101161492, "chr15": 99753195, "chr16": 96330374,
    "chr17": 84276897, "chr18": 80542538, "chr19": 61707364, "chr20": 66210255,
    "chr21": 45090682, "chr22": 51324926, "chrX": 154259566, "chrY": 62460029,
    "chrM": 16569,
}

HG38_CHROM_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    "chrM": 16569,
}


def parse_paf(paf_file, target_chr):
    """
    Parse PAF file and extract alignments for the specified chromosome.

    Returns:
        dict: {contig_name: [(target_start, target_end, strand), ...]}
        int: chromosome length from PAF
    """
    alignments = defaultdict(list)
    chrom_length = None

    with open(paf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            query_name = fields[0]
            # query_length = int(fields[1])
            # query_start = int(fields[2])
            # query_end = int(fields[3])
            strand = fields[4]
            target_name = fields[5]
            target_length = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])

            if target_name == target_chr:
                alignments[query_name].append((target_start, target_end, strand))
                if chrom_length is None:
                    chrom_length = target_length

    return alignments, chrom_length


def calculate_contig_total_length(alignments):
    """
    Calculate total mapped length for each contig.

    Returns:
        list: [(contig_name, total_length, [(start, end, strand), ...]), ...] sorted by total_length descending
    """
    contig_lengths = []
    for contig_name, aln_list in alignments.items():
        total_length = sum(end - start for start, end, _ in aln_list)
        # Sort alignments by start position
        sorted_alns = sorted(aln_list, key=lambda x: x[0])
        contig_lengths.append((contig_name, total_length, sorted_alns))

    # Sort by total length descending (longest first, will be at bottom)
    contig_lengths.sort(key=lambda x: x[1], reverse=True)
    return contig_lengths


def plot_alignments(contig_data, chrom_length, target_chr, output_file, reference_type):
    """
    Plot contig alignments.
    """
    if not contig_data:
        print(f"No alignments found for {target_chr}", file=sys.stderr)
        sys.exit(1)

    n_contigs = len(contig_data)

    # Figure setup
    fig_height = max(4, 0.4 * n_contigs + 1)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    # Colors for different strands
    color_plus = '#4477AA'   # Blue for + strand
    color_minus = '#EE6677'  # Red for - strand

    bar_height = 0.6
    total_mapped_length = 0

    # Plot from bottom (longest) to top (shortest)
    for i, (contig_name, total_length, aln_list) in enumerate(contig_data):
        y_pos = i
        total_mapped_length += total_length

        # Draw alignments
        for j, (start, end, strand) in enumerate(aln_list):
            color = color_plus if strand == '+' else color_minus
            ax.barh(y_pos, end - start, left=start, height=bar_height,
                   color=color, edgecolor='black', linewidth=0.5)

            # Draw gap markers (|) if there are gaps between alignments
            if j > 0:
                prev_end = aln_list[j-1][1]
                # Draw | at the end of previous alignment and start of current
                gap_start = prev_end
                gap_end = start

                # Only draw if there's an actual gap
                if gap_end > gap_start:
                    # Draw | markers at gap boundaries
                    marker_height = bar_height * 0.8
                    ax.plot([gap_start, gap_start],
                           [y_pos - marker_height/2, y_pos + marker_height/2],
                           color='black', linewidth=1.5)
                    ax.plot([gap_end, gap_end],
                           [y_pos - marker_height/2, y_pos + marker_height/2],
                           color='black', linewidth=1.5)

    # X-axis formatting
    ax.set_xlim(0, chrom_length)
    ax.set_xlabel(f'{target_chr} position (bp)', fontsize=10)

    # Format x-axis with Mb
    def format_mb(x, pos):
        return f'{x/1e6:.0f}Mb'
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_mb))

    # Y-axis formatting
    contig_names = [c[0] for c in contig_data]
    ax.set_yticks(range(n_contigs))
    ax.set_yticklabels(contig_names, fontsize=8)
    ax.set_ylabel('Contigs', fontsize=10)
    ax.set_ylim(-0.5, n_contigs - 0.5)

    # Title and stats
    coverage_pct = (total_mapped_length / chrom_length) * 100
    ax.set_title(f'Contig alignments to {target_chr} ({reference_type})', fontsize=12, fontweight='bold')

    # Add stats text in top right
    stats_text = f'Total length: {total_mapped_length:,} bp ({coverage_pct:.1f}%)\nn = {n_contigs} contigs'
    ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=color_plus, edgecolor='black', label='+ strand'),
        mpatches.Patch(facecolor=color_minus, edgecolor='black', label='- strand'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=8)

    # Grid
    ax.grid(axis='x', linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()

    print(f"Saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Plot contig alignments from PAF file for a specific chromosome'
    )
    parser.add_argument('paf_file', help='Input PAF file from minimap2')
    parser.add_argument('--reference', required=True, choices=['chm13', 'hg38'],
                       help='Reference genome (chm13 or hg38)')
    parser.add_argument('--chr', required=True, dest='chromosome',
                       help='Chromosome to plot (e.g., chr1, chr2, chrX)')
    parser.add_argument('-o', '--output', help='Output PNG file (default: <input>.chrN.png)')

    args = parser.parse_args()

    # Validate chromosome
    chrom_lengths = CHM13_CHROM_LENGTHS if args.reference == 'chm13' else HG38_CHROM_LENGTHS
    if args.chromosome not in chrom_lengths:
        print(f"Error: Unknown chromosome '{args.chromosome}' for {args.reference}", file=sys.stderr)
        print(f"Available: {', '.join(sorted(chrom_lengths.keys()))}", file=sys.stderr)
        sys.exit(1)

    # Parse PAF
    alignments, paf_chrom_length = parse_paf(args.paf_file, args.chromosome)

    if not alignments:
        print(f"Error: No alignments found for {args.chromosome} in {args.paf_file}", file=sys.stderr)
        sys.exit(1)

    # Use chromosome length from reference (more reliable)
    chrom_length = chrom_lengths[args.chromosome]

    # Calculate contig lengths and sort
    contig_data = calculate_contig_total_length(alignments)

    # Output file
    if args.output:
        output_file = args.output
    else:
        paf_path = Path(args.paf_file)
        output_file = str(paf_path.parent / f"{paf_path.stem}.{args.chromosome}.png")

    # Plot
    plot_alignments(contig_data, chrom_length, args.chromosome, output_file, args.reference)


if __name__ == '__main__':
    main()
