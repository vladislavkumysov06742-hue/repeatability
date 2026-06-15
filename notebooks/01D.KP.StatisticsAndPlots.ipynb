#!/usr/bin/env python3
"""
01D.KP.StatisticsAndPlots.py

Generate statistics and improved plots (histogram, heatmap, boxplot, trend, GC vs length)
from cleaned repeatability files.

Author: Vladislav Gadzhiev
Date: 12/12/2025

Usage examples:
    python 01D_statistics_plots.py --input_dir ../data/2_derived/01KP/cleaned --fasta ../data/1_raw/sequence.fasta
    python 01D_statistics_plots.py -i cleaned -f mtDNA.fasta -o figures --window 200 --max_degradation 20
    python 01D_statistics_plots.py --help
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import SeqIO
import warnings
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from scipy.stats import ttest_ind, linregress

warnings.filterwarnings('ignore')


def extract_metrics_with_gc(filepath, major_arc_start, major_arc_end, max_degradation_percent):
    """
    Extract repeatability metrics and GC content from a single file.
    Returns a dict with keys: pos, nuc, max_perfect_len, best_perfect_repeat_seq,
    best_perfect_gc, max_degraded_len.
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None

    if df.empty:
        return None

    fname = filepath.stem
    parts = fname.split('.')
    if len(parts) >= 3:
        pos = int(parts[1])
        nuc = parts[2]
    else:
        return None

    # Apply filters as in original
    mask = (
        ~((df['motif.start'] == df['repeat.start']) & (df['motif.end'] == df['repeat.end'])) &
        (df['motif.start'] >= major_arc_start) & (df['motif.end'] <= major_arc_end) &
        (df['repeat.start'] >= major_arc_start) & (df['repeat.end'] <= major_arc_end) &
        (df['motif.start'] < df['pos']) & (df['pos'] < df['motif.end'])
    )
    filtered = df[mask].copy()

    if filtered.empty:
        return {
            'pos': pos,
            'nuc': nuc,
            'max_perfect_len': np.nan,
            'best_perfect_repeat_seq': '',
            'best_perfect_gc': np.nan,
            'max_degraded_len': np.nan
        }

    # Perfect repeats (Hamming distance == 0)
    perfect = filtered[filtered['repeat.hamming.distance'] == 0]
    if not perfect.empty:
        best_perfect = perfect.loc[perfect['motif.length'].idxmax()]
        max_perfect_len = best_perfect['motif.length']
        repeat_seq = best_perfect['repeat.seq']
        gc = (repeat_seq.count('G') + repeat_seq.count('C')) / len(repeat_seq) * 100
    else:
        max_perfect_len = np.nan
        repeat_seq = ''
        gc = np.nan

    # Degraded repeats (optional)
    filtered['degradation'] = (filtered['repeat.hamming.distance'] / filtered['motif.length']) * 100
    degraded = filtered[filtered['degradation'] <= max_degradation_percent]
    max_degraded = degraded['motif.length'].max() if not degraded.empty else np.nan

    return {
        'pos': pos,
        'nuc': nuc,
        'max_perfect_len': max_perfect_len,
        'best_perfect_repeat_seq': repeat_seq,
        'best_perfect_gc': gc,
        'max_degraded_len': max_degraded
    }


def load_all_metrics(input_dir, major_arc_start, major_arc_end, max_degradation_percent, n_jobs):
    """Load metrics from all 01KP.*.txt files in input_dir using parallel processing."""
    files = list(Path(input_dir).glob("01KP.*.txt"))
    if not files:
        print(f"No files found in {input_dir}")
        return pd.DataFrame()

    print(f"Found {len(files)} files. Loading metrics with {n_jobs} workers...")
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        # Create partial function with fixed parameters
        from functools import partial
        func = partial(extract_metrics_with_gc,
                       major_arc_start=major_arc_start,
                       major_arc_end=major_arc_end,
                       max_degradation_percent=max_degradation_percent)
        metrics_list = list(executor.map(func, files))

    metrics_list = [m for m in metrics_list if m is not None]
    df = pd.DataFrame(metrics_list)
    print(f"Loaded {len(df)} records.")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Generate statistics and plots from cleaned repeatability files.",
        epilog="Example: %(prog)s -i data/2_derived/01KP/cleaned -f data/1_raw/sequence.fasta -o figures"
    )
    parser.add_argument('-i', '--input_dir', required=True,
                        help="Directory containing cleaned 01KP.*.txt files")
    parser.add_argument('-f', '--fasta', required=True,
                        help="Reference FASTA file (for reference nucleotides)")
    parser.add_argument('-o', '--output_dir', default='./figures',
                        help="Output directory for plots (default: ./figures)")
    parser.add_argument('--major_start', type=int, default=5798,
                        help="Start of major arc (default: 5798)")
    parser.add_argument('--major_end', type=int, default=16568,
                        help="End of major arc (default: 16568)")
    parser.add_argument('--max_degradation', type=int, default=20,
                        help="Maximum degradation percent for degraded metric (default: 20)")
    parser.add_argument('--window_size', type=int, default=200,
                        help="Window size for trend smoothing (default: 200)")
    parser.add_argument('--workers', type=int, default=min(mp.cpu_count(), 8),
                        help=f"Number of parallel workers (default: {min(mp.cpu_count(), 8)})")
    parser.add_argument('--dpi', type=int, default=150,
                        help="DPI for saved figures (default: 150)")
    parser.add_argument('--no_lowess', action='store_true',
                        help="Do not attempt to use LOWESS smoothing (requires statsmodels)")
    args = parser.parse_args()

    # Create output directory
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load reference sequence
    try:
        record = next(SeqIO.parse(args.fasta, "fasta"))
        mtDNA_seq = str(record.seq).upper()
        print(f"Genome length: {len(mtDNA_seq)}")
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return 1

    # Load metrics
    df_metrics = load_all_metrics(
        args.input_dir,
        args.major_start,
        args.major_end,
        args.max_degradation,
        args.workers
    )
    if df_metrics.empty:
        print("No data loaded. Exiting.")
        return 1

    # Annotate reference status
    ref_nuc = {i+1: mtDNA_seq[i] for i in range(len(mtDNA_seq))}
    df_metrics['ref_nuc'] = df_metrics['pos'].map(ref_nuc)
    df_metrics['is_ref'] = (df_metrics['nuc'] == df_metrics['ref_nuc'])

    # Prepare data for plots
    perfect = df_metrics['max_perfect_len'].dropna()
    ref_vals = df_metrics[df_metrics['is_ref']]['max_perfect_len'].dropna()
    alt_vals = df_metrics[~df_metrics['is_ref']]['max_perfect_len'].dropna()

    # Calculate change per position (max alt - ref)
    changes = []
    for pos in df_metrics['pos'].unique():
        sub = df_metrics[df_metrics['pos'] == pos]
        ref_len = sub[sub['is_ref']]['max_perfect_len'].values
        if len(ref_len) == 0 or pd.isna(ref_len[0]):
            continue
        alt_lens = sub[~sub['is_ref']]['max_perfect_len'].dropna()
        if len(alt_lens) == 0:
            continue
        max_alt = alt_lens.max()
        change = max_alt - ref_len[0]
        changes.append({'pos': pos, 'max_change': change})
    df_changes = pd.DataFrame(changes)

    # For GC plot: reference positions with non‑nan GC
    ref_with_gc = df_metrics[df_metrics['is_ref']].dropna(subset=['best_perfect_gc', 'max_perfect_len'])
    ref_with_gc = ref_with_gc[ref_with_gc['best_perfect_gc'] > 0]

    sns.set_style("whitegrid")

    # ---------- Figure 1: Histogram ----------
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.hist(perfect, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Maximum perfect repeat length (bp)')
    ax1.set_ylabel('Frequency (position‑nucleotide)')
    ax1.set_title('Distribution of maximum perfect repeat lengths')
    ax1.axvline(perfect.median(), color='red', linestyle='--', label=f'Median = {perfect.median():.1f}')
    ax1.legend()
    plt.tight_layout()
    plt.savefig(out_dir / "01_histogram_repeat_length.png", dpi=args.dpi)
    plt.close(fig1)

    # ---------- Figure 2: Heatmap ----------
    fig2, ax2 = plt.subplots(figsize=(14, 6))
    df_pivot = df_metrics.pivot_table(index='pos', columns='nuc', values='max_perfect_len', aggfunc='first')
    bin_size = 100
    bin_edges = np.arange(df_pivot.index.min(), df_pivot.index.max() + bin_size, bin_size)
    bin_labels = [f"{int(b)}-{int(b+bin_size-1)}" for b in bin_edges[:-1]]
    binned = df_pivot.groupby(pd.cut(df_pivot.index, bins=bin_edges, labels=bin_labels)).mean()
    sns.heatmap(binned.T, ax=ax2, cmap='viridis', cbar_kws={'label': 'Mean repeat length'})
    ax2.set_title('Heatmap: mean repeat length by position blocks\n(nucleotide is the one at the position)')
    ax2.set_xlabel('Position block')
    ax2.set_ylabel('Nucleotide')
    plt.tight_layout()
    plt.savefig(out_dir / "02_heatmap_repeat_by_position.png", dpi=args.dpi)
    plt.close(fig2)

    # ---------- Figure 3: Boxplot ref vs alt ----------
    fig3, ax3 = plt.subplots(figsize=(8, 6))
    data_to_plot = [ref_vals, alt_vals]
    bp = ax3.boxplot(data_to_plot, labels=['Reference nucleotide', 'Alternative nucleotides'], patch_artist=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][1].set_facecolor('lightcoral')
    ax3.set_ylabel('Maximum perfect repeat length')
    ax3.set_title('Repeatability comparison: reference vs alternative nucleotides')
    t_stat, p_val = ttest_ind(ref_vals, alt_vals, nan_policy='omit')
    ax3.text(0.5, 0.95, f't-test p = {p_val:.3e}', transform=ax3.transAxes, ha='center', fontsize=10)
    plt.tight_layout()
    plt.savefig(out_dir / "03_boxplot_ref_vs_alt.png", dpi=args.dpi)
    plt.close(fig3)

    # ---------- Figure 4: Trend of change along genome ----------
    fig4, ax4 = plt.subplots(figsize=(12, 6))
    df_changes_sorted = df_changes.sort_values('pos')
    df_changes_sorted['window'] = (df_changes_sorted['pos'] // args.window_size) * args.window_size
    window_stats = df_changes_sorted.groupby('window')['max_change'].agg(['mean', 'sem', 'count'])
    window_stats = window_stats[window_stats['count'] > 5]
    ax4.plot(window_stats.index, window_stats['mean'],
             color='darkblue', linewidth=2, label=f'Mean over {args.window_size} bp window')
    ax4.fill_between(window_stats.index,
                     window_stats['mean'] - window_stats['sem'],
                     window_stats['mean'] + window_stats['sem'],
                     color='lightblue', alpha=0.4, label='± SEM')
    ax4.axhline(y=0, color='red', linestyle='--', alpha=0.7, linewidth=1.5)
    # Optional LOWESS smoothing
    if not args.no_lowess:
        try:
            from statsmodels.nonparametric.smoothers_lowess import lowess
            lowess_result = lowess(df_changes_sorted['max_change'], df_changes_sorted['pos'], frac=0.05)
            ax4.plot(lowess_result[:, 0], lowess_result[:, 1], color='orange', linewidth=2,
                     label='LOWESS (frac=0.05)', alpha=0.8)
        except ImportError:
            print("statsmodels not available, skipping LOWESS smoothing.")
    ax4.set_xlabel('Position in mtDNA')
    ax4.set_ylabel('Mean change in repeat length (Alt - Ref)')
    ax4.set_title('Change in repeatability along the genome (window average)')
    ax4.legend(loc='upper right')
    ax4.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "04_trend_change_smoothed.png", dpi=args.dpi)
    plt.close(fig4)

    # ---------- Figure 5: GC content vs length ----------
    if not ref_with_gc.empty:
        fig5, ax5 = plt.subplots(figsize=(10, 6))
        bins_gc = np.arange(0, 101, 5)
        ref_with_gc['gc_bin'] = pd.cut(ref_with_gc['best_perfect_gc'], bins=bins_gc)
        binned_gc = ref_with_gc.groupby('gc_bin').agg(
            mean_gc=('best_perfect_gc', 'mean'),
            mean_len=('max_perfect_len', 'mean'),
            sem_len=('max_perfect_len', 'sem'),
            count=('max_perfect_len', 'count')
        ).dropna()
        binned_gc = binned_gc[binned_gc['count'] >= 5]
        if not binned_gc.empty:
            ax5.errorbar(binned_gc['mean_gc'], binned_gc['mean_len'],
                         yerr=binned_gc['sem_len'],
                         fmt='o', capsize=4, color='forestgreen',
                         ecolor='lightgreen', markersize=8, alpha=0.8,
                         label='Mean per GC bin ± SEM')
            # Linear regression on binned data
            slope, intercept, r_val, p_val, _ = linregress(binned_gc['mean_gc'], binned_gc['mean_len'])
            x_line = np.linspace(binned_gc['mean_gc'].min(), binned_gc['mean_gc'].max(), 100)
            y_line = slope * x_line + intercept
            ax5.plot(x_line, y_line, color='darkred', linewidth=2,
                     label=f'Regression: slope={slope:.2f}, R²={r_val**2:.3f}, p={p_val:.2e}')
            ax5.set_xlabel('Repeat GC content (%)')
            ax5.set_ylabel('Mean repeat length (bp)')
            ax5.set_title('Dependence of repeat length on its GC content (binned)')
            ax5.legend()
            ax5.grid(True, linestyle='--', alpha=0.3)
            plt.tight_layout()
            plt.savefig(out_dir / "05_binned_gc_vs_length.png", dpi=args.dpi)
            plt.close(fig5)
        else:
            print("Not enough data for GC vs length plot (bins with ≥5 observations).")
    else:
        print("No reference positions with positive GC content found, skipping Figure 5.")

    print(f"\n✅ All plots saved to {out_dir}")
    print("   - 01_histogram_repeat_length.png")
    print("   - 02_heatmap_repeat_by_position.png")
    print("   - 03_boxplot_ref_vs_alt.png")
    print("   - 04_trend_change_smoothed.png")
    if not ref_with_gc.empty and 'binned_gc' in locals() and not binned_gc.empty:
        print("   - 05_binned_gc_vs_length.png")
    return 0


if __name__ == "__main__":
    exit(main())
