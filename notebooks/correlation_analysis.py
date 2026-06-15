#!/usr/bin/env python3
"""
correlation_analysis.py

Compare repeatability profiles between two species (e.g., human and pika) using
cleaned repeatability files. Computes overall Pearson correlation and sliding‑window
correlation, then generates a three‑panel figure and saves results.

Author: Vladislav Gadziev
Date: 13/12/2025

Usage examples:
    python correlation_analysis.py --human_dir data/human/cleaned --pika_dir data/pika/cleaned
    python correlation_analysis.py -H human/cleaned -P pika/cleaned -o comparison --window 200 --dpi 150
    python correlation_analysis.py --help
"""

import argparse
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')


def get_genome_length_from_files(results_dir):
    """
    Determine genome length as the maximum position number found in any
    01KP.*.txt file in results_dir.
    """
    max_pos = 0
    pattern = re.compile(r'01KP\.(\d+)\.([ATCG])\.txt$')
    for fpath in Path(results_dir).glob("01KP.*.txt"):
        match = pattern.match(fpath.name)
        if match:
            pos = int(match.group(1))
            if pos > max_pos:
                max_pos = pos
    if max_pos == 0:
        raise ValueError(f"No 01KP.*.txt files found in {results_dir}")
    return max_pos


def load_repeatability_for_genome(results_dir, genome_len=None, nucleotides=None):
    """
    Load all 01KP.*.txt files from results_dir and return a 1D array (length = genome_len)
    containing the maximum repeatability (effective length of perfect repeats) for each position.

    For each position, we take the maximum across all four nucleotides.
    """
    if nucleotides is None:
        nucleotides = ['A', 'T', 'G', 'C']
    if genome_len is None:
        genome_len = get_genome_length_from_files(results_dir)

    nuc_to_idx = {nuc: i for i, nuc in enumerate(nucleotides)}
    # matrix: positions x nucleotides
    repeatability = np.zeros((genome_len, len(nucleotides)), dtype=float)

    files = list(Path(results_dir).glob("01KP.*.txt"))
    for fpath in tqdm(files, desc=f"Loading {Path(results_dir).name}"):
        parts = fpath.stem.split('.')
        if len(parts) != 3 or parts[0] != "01KP":
            continue
        try:
            pos = int(parts[1])
            nuc = parts[2]
            if nuc not in nuc_to_idx:
                continue
        except Exception:
            continue
        if pos < 1 or pos > genome_len:
            continue
        try:
            df = pd.read_csv(fpath, sep='\t')
        except Exception:
            continue
        # Use perfect repeats (hamming distance == 0)
        perfect = df[df['repeat.hamming.distance'] == 0]
        if not perfect.empty:
            best = perfect['effective.length'].max()
            repeatability[pos-1, nuc_to_idx[nuc]] = best
        # otherwise remains 0

    # For each position, take the maximum across nucleotides
    profile = repeatability.max(axis=1)
    return profile


def rolling_correlation(x, y, window):
    """Compute sliding‑window Pearson correlation between two arrays."""
    if len(x) != len(y):
        raise ValueError("Arrays must have the same length")
    n = len(x)
    corrs = np.full(n, np.nan)
    half = window // 2
    for i in range(half, n - half):
        start = i - half
        end = i + half + 1
        x_win = x[start:end]
        y_win = y[start:end]
        mask = ~(np.isnan(x_win) | np.isnan(y_win))
        if np.sum(mask) >= 3:
            corr, _ = pearsonr(x_win[mask], y_win[mask])
            corrs[i] = corr
    return corrs


def main():
    parser = argparse.ArgumentParser(
        description="Compare repeatability profiles between two species (e.g., human and pika).",
        epilog="Example: %(prog)s --human_dir data/human/cleaned --pika_dir data/pika/cleaned --output comparison"
    )
    parser.add_argument('-H', '--human_dir', required=True,
                        help="Directory containing cleaned 01KP.*.txt files for human (or first genome)")
    parser.add_argument('-P', '--pika_dir', required=True,
                        help="Directory containing cleaned 01KP.*.txt files for pika (or second genome)")
    parser.add_argument('-o', '--output_dir', default='./comparison_results',
                        help="Output directory for figures and CSV (default: ./comparison_results)")
    parser.add_argument('--window', type=int, default=100,
                        help="Window size for sliding correlation (in nucleotides, default: 100)")
    parser.add_argument('--dpi', type=int, default=150,
                        help="DPI for output figure (default: 150)")
    parser.add_argument('--nucleotides', nargs='+', default=['A','T','G','C'],
                        help="Nucleotides to consider (default: A T G C)")
    parser.add_argument('--skip_plot', action='store_true',
                        help="Skip generating the figure (only save CSV)")
    args = parser.parse_args()

    # Create output directory
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Determine genome lengths (common length = min of the two)
    print("Determining genome lengths from files...")
    human_len = get_genome_length_from_files(args.human_dir)
    pika_len = get_genome_length_from_files(args.pika_dir)
    print(f"Human genome length (from files): {human_len} bp")
    print(f"Pika genome length (from files): {pika_len} bp")
    common_len = min(human_len, pika_len)
    print(f"Using common length for comparison: {common_len} bp")

    # Load profiles
    print("Loading human repeatability profile...")
    human_profile = load_repeatability_for_genome(args.human_dir, genome_len=common_len,
                                                  nucleotides=args.nucleotides)
    print("Loading pika repeatability profile...")
    pika_profile = load_repeatability_for_genome(args.pika_dir, genome_len=common_len,
                                                 nucleotides=args.nucleotides)

    # Overall correlation (only positions where both > 0)
    mask = (human_profile > 0) & (pika_profile > 0)
    if np.sum(mask) > 1:
        overall_corr, p_value = pearsonr(human_profile[mask], pika_profile[mask])
        print(f"\nOverall Pearson correlation (positions with repeats in both genomes): "
              f"r = {overall_corr:.4f}, p = {p_value:.2e}")
    else:
        overall_corr = np.nan
        p_value = np.nan
        print("Insufficient data for overall correlation.")

    # Sliding window correlation
    print(f"Computing sliding correlation (window = {args.window} bp)...")
    rolling_corr = rolling_correlation(human_profile, pika_profile, args.window)

    # Save summary CSV
    summary_df = pd.DataFrame({
        'position': np.arange(1, common_len + 1),
        'human_repeatability': human_profile,
        'pika_repeatability': pika_profile,
        'rolling_correlation': rolling_corr
    })
    csv_path = out_dir / 'repeatability_comparison.csv'
    summary_df.to_csv(csv_path, index=False)
    print(f"Saved summary CSV to {csv_path}")

    # Generate figure (three panels)
    if not args.skip_plot:
        fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True, dpi=args.dpi)

        # Panel 1: repeatability profiles
        axes[0].plot(human_profile, label='Human', alpha=0.7, color='blue')
        axes[0].plot(pika_profile, label='Pika', alpha=0.7, color='orange')
        axes[0].set_ylabel('Max repeatability\n(effective length)')
        axes[0].legend()
        axes[0].set_title('Repeatability profiles (max over nucleotides)')
        axes[0].grid(alpha=0.3)

        # Panel 2: rolling correlation
        axes[1].plot(rolling_corr, color='green')
        axes[1].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        axes[1].set_ylabel(f'Rolling correlation (window={args.window} bp)')
        axes[1].set_ylim(-1, 1)
        axes[1].grid(alpha=0.3)
        axes[1].set_title(f'Sliding Pearson correlation (window {args.window} bp)')

        # Panel 3: scatter plot
        axes[2].scatter(human_profile, pika_profile, alpha=0.3, s=1)
        if not np.isnan(overall_corr):
            # Fit linear regression line
            z = np.polyfit(human_profile, pika_profile, 1)
            p_line = np.poly1d(z)
            x_max = max(human_profile.max(), pika_profile.max())
            x_line = np.linspace(0, x_max, 100)
            axes[2].plot(x_line, p_line(x_line), "r--", label=f'r = {overall_corr:.3f}')
        axes[2].set_xlabel('Human repeatability')
        axes[2].set_ylabel('Pika repeatability')
        axes[2].legend()
        axes[2].grid(alpha=0.3)
        axes[2].set_title('Position‑wise repeatability comparison')
        axes[2].set_xlim(0, None)
        axes[2].set_ylim(0, None)

        plt.tight_layout()
        fig_path = out_dir / 'repeatability_comparison.png'
        plt.savefig(fig_path, dpi=args.dpi, bbox_inches='tight')
        plt.close()
        print(f"Figure saved to {fig_path}")
    else:
        print("Plot generation skipped.")

    print(f"\nAll results saved in {out_dir}")
    return 0


if __name__ == "__main__":
    exit(main())
