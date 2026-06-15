#!/usr/bin/env python3
"""
scatter_plot.py

Generate scatter plot (reference repeat length vs |ΔR|) and histogram of ΔR
for all positions or user‑specified positions. Also saves ΔR values to CSV.

Author: Vladislav Gadziev
Date: 16/12/2025

Usage examples:
    python scatter_plot.py --input_dir data/2_derived/01KP/cleaned --fasta data/1_raw/sequence.fasta
    python scatter_plot.py -i cleaned -f mtDNA.fasta -o my_figure.png --positions 8251 8473 12705 --dpi 300
    python scatter_plot.py --help
"""

import argparse
import sys
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from Bio import SeqIO

warnings.filterwarnings('ignore')


# ----------------------------------------------------------------------
# Helper functions (copied from original notebook)
# ----------------------------------------------------------------------
def is_transition(ref, alt):
    transitions = {('A','G'), ('G','A'), ('C','T'), ('T','C')}
    return (ref, alt) in transitions


def load_and_aggregate(data_dir, ref_nuc_dict):
    """
    Read all 01KP.*.txt files, compute per‑allele summary statistics.
    Returns DataFrame with columns:
        position, nucleotide, is_reference, mean_length, max_length, median_length,
        mean_effective, max_effective, mean_gc, median_gc, perfect_repeats,
        total_repeats, percent_perfect, total_bases, avg_sequence_length, file_path
    """
    results = []
    data_path = Path(data_dir)
    pattern = re.compile(r"01KP\.(\d+)\.([ATGC])\.txt")
    files = list(data_path.glob("01KP.*.txt"))
    print(f"Found {len(files)} files.")

    for fpath in files:
        match = pattern.match(fpath.name)
        if not match:
            continue
        pos = int(match.group(1))
        nuc = match.group(2)

        try:
            df = pd.read_csv(fpath, sep='\t')
            if df.empty:
                continue

            mean_length = df['motif.length'].mean()
            max_length = df['motif.length'].max()
            median_length = df['motif.length'].median()

            if 'effective.length' in df.columns:
                mean_effective = df['effective.length'].mean()
                max_effective = df['effective.length'].max()
            else:
                mean_effective = mean_length
                max_effective = max_length

            gc_values = []
            for seq in df['motif.seq']:
                if isinstance(seq, str):
                    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
                    gc_values.append(gc)
            mean_gc = np.mean(gc_values) if gc_values else np.nan
            median_gc = np.median(gc_values) if gc_values else np.nan

            perfect = len(df[df['repeat.hamming.distance'] == 0]) if 'repeat.hamming.distance' in df.columns else np.nan
            total = len(df)
            percent_perfect = (perfect / total * 100) if total > 0 else 0

            total_bases = df['motif.length'].sum() if total > 0 else 0
            avg_seq_len = total_bases / total if total > 0 else 0

            results.append({
                'position': pos,
                'nucleotide': nuc,
                'is_reference': (ref_nuc_dict.get(pos) == nuc),
                'mean_length': mean_length,
                'max_length': max_length,
                'median_length': median_length,
                'mean_effective': mean_effective,
                'max_effective': max_effective,
                'mean_gc': mean_gc,
                'median_gc': median_gc,
                'perfect_repeats': perfect,
                'total_repeats': total,
                'percent_perfect': percent_perfect,
                'total_bases': total_bases,
                'avg_sequence_length': avg_seq_len,
                'file_path': str(fpath)
            })
        except Exception as e:
            print(f"Error reading {fpath.name}: {e}")

    return pd.DataFrame(results)


def calculate_delta_R(df_all, ref_nuc_dict, positions=None):
    """
    Compute ΔR for every alternative allele relative to the reference.
    If positions is given, restrict to those positions.
    Returns DataFrame with columns: position, ref, alt, ref_length, alt_length,
    delta_R, abs_delta_R, delta_R_effective, delta_GC, delta_perfect,
    mutation, is_transition, ref_perfect, alt_perfect.
    """
    if positions is not None:
        df_all = df_all[df_all['position'].isin(positions)]
    delta_data = []
    for pos in df_all['position'].unique():
        pos_data = df_all[df_all['position'] == pos]
        ref_nuc = ref_nuc_dict.get(pos)
        if not ref_nuc:
            continue
        ref_row = pos_data[pos_data['nucleotide'] == ref_nuc]
        if ref_row.empty:
            continue
        ref_row = ref_row.iloc[0]
        for _, row in pos_data.iterrows():
            if row['nucleotide'] != ref_nuc:
                delta_data.append({
                    'position': pos,
                    'ref': ref_nuc,
                    'alt': row['nucleotide'],
                    'ref_length': ref_row['mean_length'],
                    'alt_length': row['mean_length'],
                    'delta_R': row['mean_length'] - ref_row['mean_length'],
                    'abs_delta_R': abs(row['mean_length'] - ref_row['mean_length']),
                    'delta_R_effective': row['mean_effective'] - ref_row['mean_effective'],
                    'delta_GC': row['mean_gc'] - ref_row['mean_gc'],
                    'delta_perfect': row['percent_perfect'] - ref_row['percent_perfect'],
                    'mutation': f"{ref_nuc}→{row['nucleotide']}",
                    'is_transition': is_transition(ref_nuc, row['nucleotide']),
                    'ref_perfect': ref_row['percent_perfect'],
                    'alt_perfect': row['percent_perfect']
                })
    return pd.DataFrame(delta_data)


# ----------------------------------------------------------------------
# Plotting function
# ----------------------------------------------------------------------
def plot_scatter_histogram(df_delta, output_path, dpi=300):
    """Generate scatter plot (ref_length vs |ΔR|, colored by ΔGC) and histogram of ΔR."""
    if df_delta.empty:
        print("No data to plot.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=dpi)

    # Scatter plot: reference length vs |ΔR|
    sc = axes[0].scatter(df_delta['ref_length'], df_delta['abs_delta_R'],
                         c=df_delta['delta_GC'], cmap='coolwarm',
                         alpha=0.7, s=80, edgecolor='black', linewidth=0.5)
    axes[0].set_xlabel('Mean reference repeat length (bp)', fontsize=11)
    axes[0].set_ylabel('|ΔR|', fontsize=11)
    axes[0].set_title('Dependence of |ΔR| on reference repeat length', fontsize=13, fontweight='bold')
    axes[0].grid(True, alpha=0.3)

    # Linear regression if enough variation
    if df_delta['ref_length'].nunique() > 1 and df_delta['abs_delta_R'].nunique() > 1:
        slope, intercept, r_val, p_val, _ = stats.linregress(df_delta['ref_length'], df_delta['abs_delta_R'])
        x_line = np.array([df_delta['ref_length'].min(), df_delta['ref_length'].max()])
        y_line = intercept + slope * x_line
        axes[0].plot(x_line, y_line, 'r--', linewidth=2, label=f'R²={r_val**2:.3f}, p={p_val:.3e}')
        axes[0].legend()

    cbar = plt.colorbar(sc, ax=axes[0])
    cbar.set_label('ΔGC (%)', fontsize=10)

    # Histogram of ΔR
    axes[1].hist(df_delta['delta_R'], bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    axes[1].axvline(x=0, color='red', linestyle='--', linewidth=2, alpha=0.7)
    axes[1].set_xlabel('ΔR', fontsize=11)
    axes[1].set_ylabel('Frequency', fontsize=11)
    axes[1].set_title('Distribution of ΔR', fontsize=13, fontweight='bold')
    axes[1].grid(True, alpha=0.3)

    mean_delta = df_delta['delta_R'].mean()
    std_delta = df_delta['delta_R'].std()
    axes[1].text(0.95, 0.95, f'μ = {mean_delta:.3f}\nσ = {std_delta:.3f}\nn = {len(df_delta)}',
                 transform=axes[1].transAxes, fontsize=10,
                 verticalalignment='top', horizontalalignment='right',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.suptitle('Effect of nucleotide substitutions on repeat length', fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Figure saved: {output_path}")


# ----------------------------------------------------------------------
# Main CLI
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate scatter plot (ref length vs |ΔR|) and histogram of ΔR.",
        epilog="Example: %(prog)s --input_dir data/cleaned --fasta mtDNA.fasta --output my_plot.png"
    )
    parser.add_argument('-i', '--input_dir', required=True,
                        help="Directory containing cleaned 01KP.*.txt files")
    parser.add_argument('-f', '--fasta', required=True,
                        help="Reference FASTA file (used to determine reference nucleotides)")
    parser.add_argument('-o', '--output', default='scatter_histogram.png',
                        help="Output figure file name (default: scatter_histogram.png)")
    parser.add_argument('-p', '--positions', type=int, nargs='+',
                        help="List of positions to analyze (if not given, use all positions with files)")
    parser.add_argument('--csv', default='delta_R_values.csv',
                        help="Output CSV file for ΔR values (default: delta_R_values.csv)")
    parser.add_argument('--dpi', type=int, default=300,
                        help="DPI for output figure (default: 300)")
    args = parser.parse_args()

    # Load reference genome and build ref_nuc dictionary
    try:
        record = next(SeqIO.parse(args.fasta, 'fasta'))
        ref_seq = str(record.seq).upper()
        ref_nuc_dict = {i+1: ref_seq[i] for i in range(len(ref_seq))}
        print(f"Reference genome length: {len(ref_seq)}")
    except Exception as e:
        print(f"Error reading FASTA: {e}")
        return 1

    # Load and aggregate data
    print("Loading repeatability files...")
    df_all = load_and_aggregate(args.input_dir, ref_nuc_dict)
    if df_all.empty:
        print("No data loaded. Exiting.")
        return 1
    print(f"Loaded {len(df_all)} allele records from {df_all['position'].nunique()} positions.")

    # Compute ΔR (optionally subset positions)
    if args.positions:
        positions_set = set(args.positions)
        available = set(df_all['position'].unique())
        missing = positions_set - available
        if missing:
            print(f"Warning: positions {sorted(missing)} not found in data. Skipping them.")
        positions_to_use = [p for p in args.positions if p in available]
        if not positions_to_use:
            print("No valid positions specified. Exiting.")
            return 1
        print(f"Analyzing {len(positions_to_use)} specified positions.")
    else:
        positions_to_use = None
        print("Analyzing all positions.")

    df_delta = calculate_delta_R(df_all, ref_nuc_dict, positions=positions_to_use)
    if df_delta.empty:
        print("No ΔR values computed. Exiting.")
        return 1
    print(f"Computed ΔR for {len(df_delta)} mutations.")

    # Save ΔR to CSV
    csv_path = Path(args.csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df_delta.to_csv(csv_path, index=False)
    print(f"ΔR values saved to {csv_path}")

    # Generate plot
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_scatter_histogram(df_delta, output_path, dpi=args.dpi)

    # Print summary
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)
    print(f"Total positions analyzed: {df_delta['position'].nunique()}")
    print(f"Total mutations: {len(df_delta)}")
    print(f"Mean ΔR: {df_delta['delta_R'].mean():.4f}  median: {df_delta['delta_R'].median():.4f}")
    inc = (df_delta['delta_R'] > 0).sum()
    dec = (df_delta['delta_R'] < 0).sum()
    print(f"Increased: {inc} ({inc/len(df_delta)*100:.1f}%)  Decreased: {dec} ({dec/len(df_delta)*100:.1f}%)")
    return 0


if __name__ == "__main__":
    exit(main())
