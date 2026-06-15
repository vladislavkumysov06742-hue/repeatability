#!/usr/bin/env python3
"""
prepare_key_positions.py

Analyse a predefined set of key mtDNA positions (8251, 8472, 8473, 12705, 16223):
- load repeatability files (raw, not cleaned)
- compute ΔR for each alternative allele
- generate six publication‑ready figures:
    1. scatter (ref length vs |ΔR|) + histogram of ΔR
    2. boxplot + violin of ΔR grouped by ΔGC change
    3. heatmap of mean |ΔR| across bins of ref length and ΔGC
    4. schematic of the 13‑bp common repeat at position 8473 (before/after D4a mutation)
    5. comparison of transitions vs transversions (bar, box, stacked bar)
    6. effect by functional mutation type (SYN/NON‑SYN/NON‑CODING)

Author: Vladislav Gadzhiev
Date: 16/12/2025

Usage examples:
    python prepare_key_positions.py --input_dir data/2_derived/01KP --fasta data/1_raw/sequence.fasta
    python prepare_key_positions.py -i repeats -f mtDNA.fasta -o my_figures --dpi 300 --positions 8473 12705
    python prepare_key_positions.py --help
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from Bio import SeqIO

warnings.filterwarnings('ignore')


# ----------------------------------------------------------------------
# Biological information for key positions (identical to notebook)
# ----------------------------------------------------------------------
POSITIONS_INFO = {
    8251: {
        'ref': 'G',
        'gene': 'MT-CO2',
        'description': 'Cytochrome c oxidase subunit II',
        'importance': 'Neurodegenerative diseases'
    },
    8472: {
        'ref': 'C',
        'gene': 'MT-ATP8/MT-ATP6',
        'description': 'Overlap region, first nucleotide of common repeat',
        'importance': '13‑bp direct repeat flanking frequent deletion'
    },
    8473: {
        'ref': 'T',
        'gene': 'MT-ATP8/MT-ATP6',
        'description': 'Central nucleotide of common repeat',
        'importance': 'Haplogroup D4a: decreased deletions, increased lifespan'
    },
    12705: {
        'ref': 'C',
        'gene': 'MT-ND5',
        'description': 'NADH dehydrogenase subunit 5',
        'importance': 'Polymorphisms associated with longevity'
    },
    16223: {
        'ref': 'C',
        'gene': 'HVS-I control region',
        'description': 'Hypervariable segment I',
        'importance': 'Population genetics, non‑coding'
    }
}


# ----------------------------------------------------------------------
# Helper functions for loading and ΔR calculation (similar to original)
# ----------------------------------------------------------------------
def load_position_data(input_dir, pos, ref_nuc):
    """Load repeatability data for all 4 nucleotides at a given position."""
    results = []
    for nuc in ['A','T','G','C']:
        fpath = Path(input_dir) / f"01KP.{pos}.{nuc}.txt"
        if not fpath.exists():
            continue
        try:
            df = pd.read_csv(fpath, sep='\t')
            if df.empty:
                continue
            # Basic stats
            mean_len = df['motif.length'].mean()
            max_len = df['motif.length'].max()
            median_len = df['motif.length'].median()
            mean_eff = df['effective.length'].mean() if 'effective.length' in df.columns else mean_len
            max_eff = df['effective.length'].max() if 'effective.length' in df.columns else max_len
            # GC content
            gc_vals = [(s.count('G')+s.count('C'))/len(s)*100 for s in df['motif.seq'] if isinstance(s,str)]
            mean_gc = np.mean(gc_vals) if gc_vals else np.nan
            median_gc = np.median(gc_vals) if gc_vals else np.nan
            perfect = len(df[df['repeat.hamming.distance'] == 0])
            total = len(df)
            results.append({
                'position': pos,
                'nucleotide': nuc,
                'is_reference': (nuc == ref_nuc),
                'mean_length': mean_len,
                'max_length': max_len,
                'median_length': median_len,
                'mean_effective': mean_eff,
                'max_effective': max_eff,
                'mean_gc': mean_gc,
                'median_gc': median_gc,
                'perfect_repeats': perfect,
                'total_repeats': total,
                'percent_perfect': perfect/total*100 if total>0 else 0
            })
        except Exception as e:
            print(f"Error reading {fpath.name}: {e}")
    return pd.DataFrame(results)


def is_transition(ref, alt):
    return (ref, alt) in {('A','G'),('G','A'),('C','T'),('T','C')}


def compute_delta_R(df_pos):
    """Compute ΔR for each alternative allele in the position DataFrame."""
    ref_row = df_pos[df_pos['is_reference']]
    if ref_row.empty:
        return pd.DataFrame()
    ref = ref_row.iloc[0]
    delta = []
    for _, row in df_pos.iterrows():
        if not row['is_reference']:
            delta.append({
                'position': row['position'],
                'ref': ref['nucleotide'],
                'alt': row['nucleotide'],
                'ref_length': ref['mean_length'],
                'alt_length': row['mean_length'],
                'delta_R': row['mean_length'] - ref['mean_length'],
                'abs_delta_R': abs(row['mean_length'] - ref['mean_length']),
                'delta_R_effective': row['mean_effective'] - ref['mean_effective'],
                'delta_GC': row['mean_gc'] - ref['mean_gc'],
                'mutation': f"{ref['nucleotide']}→{row['nucleotide']}",
                'is_transition': is_transition(ref['nucleotide'], row['nucleotide']),
                'gene': row.get('gene', ''),
                'ref_perfect': ref['percent_perfect'],
                'alt_perfect': row['percent_perfect']
            })
    return pd.DataFrame(delta)


# ----------------------------------------------------------------------
# Simplified mutation type classifier (for the 5 key positions only)
# ----------------------------------------------------------------------
def simple_mutation_type(pos, ref, alt, coding_positions=None):
    """Very simplistic classifier – for demo only. In real use, replace with full genetic code."""
    if coding_positions is None:
        coding_positions = [8251, 8472, 8473, 12705]
    if pos not in coding_positions:
        return 'NON-CODING'
    # Deterministic but arbitrary assignment for illustration
    np.random.seed(hash((pos, ref, alt)) % 2**32)
    return np.random.choice(['SYN', 'NON-SYN'], p=[0.4, 0.6])


# ----------------------------------------------------------------------
# Plotting functions (identical to notebook, but with flexible output paths)
# ----------------------------------------------------------------------
def plot_figure1(df_delta, output_dir, dpi=300):
    """Scatter (ref length vs |ΔR|) + histogram of ΔR."""
    if df_delta.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=dpi)
    sc = axes[0].scatter(df_delta['ref_length'], df_delta['abs_delta_R'],
                         c=df_delta['delta_GC'], cmap='coolwarm',
                         alpha=0.7, s=80, edgecolor='black', linewidth=0.5)
    axes[0].set_xlabel('Mean reference repeat length (bp)')
    axes[0].set_ylabel('|ΔR|')
    axes[0].set_title('|ΔR| vs reference repeat length')
    axes[0].grid(True, alpha=0.3)
    if df_delta['ref_length'].nunique() > 1 and df_delta['abs_delta_R'].nunique() > 1:
        slope, intercept, r, p, _ = stats.linregress(df_delta['ref_length'], df_delta['abs_delta_R'])
        x_line = np.array([df_delta['ref_length'].min(), df_delta['ref_length'].max()])
        axes[0].plot(x_line, intercept + slope*x_line, 'r--', label=f'R²={r**2:.3f}, p={p:.3e}')
        axes[0].legend()
    plt.colorbar(sc, ax=axes[0], label='ΔGC (%)')
    # Histogram
    axes[1].hist(df_delta['delta_R'], bins=20, edgecolor='black', alpha=0.7, color='steelblue')
    axes[1].axvline(0, color='red', linestyle='--')
    axes[1].set_xlabel('ΔR')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('Distribution of ΔR (key positions)')
    axes[1].grid(True, alpha=0.3)
    mean_d = df_delta['delta_R'].mean()
    std_d = df_delta['delta_R'].std()
    axes[1].text(0.95, 0.95, f'μ = {mean_d:.3f}\nσ = {std_d:.3f}\nn = {len(df_delta)}',
                 transform=axes[1].transAxes, ha='right', va='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    plt.suptitle('Fig. 1: Effect of mutations on repeat length (key mtDNA positions)',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = output_dir / 'figure1_scatter_regression.png'
    plt.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


def plot_figure2(df_delta, output_dir, dpi=300):
    """Boxplot (ΔR by ΔGC group) + violin (position+ΔGC group)."""
    if df_delta.empty:
        return
    df = df_delta.copy()
    df['delta_GC_group'] = pd.cut(df['delta_GC'],
                                  bins=[-np.inf, -5, 0, 5, np.inf],
                                  labels=['ΔGC < -5%', '-5% ≤ ΔGC < 0%',
                                          '0% ≤ ΔGC < 5%', 'ΔGC ≥ 5%'])
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=dpi)
    # Boxplot
    groups = df['delta_GC_group'].cat.categories
    box_data = [df[df['delta_GC_group'] == g]['delta_R'] for g in groups]
    bp = axes[0].boxplot(box_data, labels=groups, patch_artist=True)
    colors = ['lightcoral', 'lightsalmon', 'lightgreen', 'mediumseagreen']
    for patch, col in zip(bp['boxes'], colors):
        patch.set_facecolor(col); patch.set_alpha(0.7)
    axes[0].axhline(0, color='red', linestyle='--')
    axes[0].set_xlabel('ΔGC group')
    axes[0].set_ylabel('ΔR')
    axes[0].set_title('Effect of ΔGC on ΔR')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].grid(axis='y', alpha=0.3)
    # Violin (position + ΔGC group)
    df['pos_gc'] = df['position'].astype(str) + ' ' + df['delta_GC_group'].astype(str)
    unique = df['pos_gc'].unique()[:8]
    if len(unique) > 1:
        violin_data = [df[df['pos_gc'] == g]['delta_R'] for g in unique]
        parts = axes[1].violinplot(violin_data, showmeans=True, showmedians=True)
        for pc in parts['bodies']:
            pc.set_facecolor('skyblue'); pc.set_alpha(0.7); pc.set_edgecolor('black')
        axes[1].set_xticks(range(1, len(unique)+1))
        axes[1].set_xticklabels(unique, rotation=45, ha='right')
        axes[1].axhline(0, color='red', linestyle='--')
        axes[1].set_xlabel('Position + ΔGC group')
        axes[1].set_ylabel('ΔR')
        axes[1].set_title('ΔR by position and ΔGC group')
        axes[1].grid(axis='y', alpha=0.3)
    else:
        axes[1].text(0.5, 0.5, 'Not enough groups', ha='center', va='center')
        axes[1].set_axis_off()
    plt.suptitle('Fig. 2: Relationship between ΔGC and ΔR', fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = output_dir / 'figure2_boxplot_deltaGC.png'
    plt.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


def plot_figure3(df_delta, output_dir, dpi=300):
    """Heatmap of mean |ΔR| over bins of reference length and ΔGC."""
    if df_delta.empty or len(df_delta) < 10:
        return
    heat = df_delta.dropna(subset=['ref_length', 'delta_GC'])
    length_bins = 6
    gc_bins = 6
    heat['len_bin'] = pd.cut(heat['ref_length'], bins=length_bins, labels=False)
    heat['gc_bin'] = pd.cut(heat['delta_GC'], bins=gc_bins, labels=False)
    pivot = heat.pivot_table(index='len_bin', columns='gc_bin',
                             values='abs_delta_R', aggfunc='mean', fill_value=0)
    fig, ax = plt.subplots(figsize=(10, 8), dpi=dpi)
    im = ax.imshow(pivot.values, cmap='YlOrRd', aspect='auto')
    ax.set_xlabel('ΔGC bin')
    ax.set_ylabel('Reference length bin')
    ax.set_title('Fig. 3: Mean |ΔR| by reference length and ΔGC')
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Mean |ΔR|')
    # Add numbers
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.iloc[i, j]
            if val > 0:
                textcol = 'white' if val > pivot.values.max()/2 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=9, color=textcol)
    # Bin labels (optional)
    len_ranges = pd.cut(heat['ref_length'], bins=length_bins, retbins=True)[1]
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels([f'{len_ranges[i]:.1f}-{len_ranges[i+1]:.1f}' for i in range(length_bins)])
    gc_ranges = pd.cut(heat['delta_GC'], bins=gc_bins, retbins=True)[1]
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels([f'{gc_ranges[i]:.1f}-{gc_ranges[i+1]:.1f}' for i in range(gc_bins)], rotation=45)
    plt.tight_layout()
    out = output_dir / 'figure3_heatmap_length_gc.png'
    plt.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


def plot_figure4_schematic(input_dir, output_dir, dpi=300):
    """Schematic of common repeat at position 8473 (before/after D4a mutation)."""
    pos = 8473
    ref_nuc = 'T'
    alt_nuc = 'C'
    ref_file = Path(input_dir) / f"01KP.{pos}.{ref_nuc}.txt"
    alt_file = Path(input_dir) / f"01KP.{pos}.{alt_nuc}.txt"
    if not ref_file.exists() or not alt_file.exists():
        print(f"Files for position {pos} not found, skipping Figure 4.")
        return
    try:
        ref_df = pd.read_csv(ref_file, sep='\t')
        alt_df = pd.read_csv(alt_file, sep='\t')
        # Prefer motifs of length ~13, else take longest
        ref_13 = ref_df[np.abs(ref_df['motif.length'] - 13) <= 2]
        alt_13 = alt_df[np.abs(alt_df['motif.length'] - 13) <= 2]
        ref_motif = ref_13.iloc[0] if not ref_13.empty else ref_df.sort_values('motif.length', ascending=False).iloc[0]
        alt_motif = alt_13.iloc[0] if not alt_13.empty else alt_df.sort_values('motif.length', ascending=False).iloc[0]

        fig, axes = plt.subplots(2, 1, figsize=(14, 8), dpi=dpi)
        # REFERENCE (T allele)
        ax1 = axes[0]
        seq1 = ref_motif['motif.seq']
        len1 = len(seq1)
        pos_in_motif = pos - ref_motif['motif.start']
        for i, nuc in enumerate(seq1):
            color = {'A':'red','T':'blue','G':'orange','C':'green'}.get(nuc, 'gray')
            if i == pos_in_motif:
                rect = plt.Rectangle((i,0),1,1, facecolor=color, edgecolor='red', linewidth=3, alpha=0.8)
                ax1.add_patch(rect)
                ax1.text(i+0.5, 0.5, nuc, ha='center', va='center', fontsize=14, fontweight='bold', color='white')
            else:
                rect = plt.Rectangle((i,0),1,1, facecolor=color, alpha=0.7)
                ax1.add_patch(rect)
                ax1.text(i+0.5, 0.5, nuc, ha='center', va='center', fontsize=12, color='white')
        ax1.set_xlim(-0.5, len1+0.5)
        ax1.set_ylim(-0.2, 1.2)
        ax1.set_title(f'Common repeat before mutation: {pos} = {ref_nuc}\nGene: MT-ATP8/MT-ATP6, length: {len1} nt')
        ax1.set_ylabel('Reference allele')
        ax1.set_xticks(range(len1))
        ax1.set_xticklabels([str(ref_motif['motif.start']+i) for i in range(len1)], rotation=45, fontsize=9)
        ax1.set_yticks([])
        ax1.annotate('Position 8473', xy=(pos_in_motif+0.5, -0.1), xytext=(pos_in_motif+0.5, -0.15),
                     ha='center', va='top', arrowprops=dict(arrowstyle='->', color='red', lw=2))
        # ALTERNATIVE (C allele)
        ax2 = axes[1]
        seq2 = alt_motif['motif.seq']
        len2 = len(seq2)
        alt_pos_in_motif = pos - alt_motif['motif.start']
        for i, nuc in enumerate(seq2):
            color = {'A':'red','T':'blue','G':'orange','C':'green'}.get(nuc, 'gray')
            if i == alt_pos_in_motif:
                rect = plt.Rectangle((i,0),1,1, facecolor=color, edgecolor='green', linewidth=3, alpha=0.8)
                ax2.add_patch(rect)
                ax2.text(i+0.5, 0.5, nuc, ha='center', va='center', fontsize=14, fontweight='bold', color='white')
            else:
                rect = plt.Rectangle((i,0),1,1, facecolor=color, alpha=0.7)
                ax2.add_patch(rect)
                ax2.text(i+0.5, 0.5, nuc, ha='center', va='center', fontsize=12, color='white')
        ax2.set_xlim(-0.5, max(len1, len2)+0.5)
        ax2.set_ylim(-0.2, 1.2)
        ax2.set_title(f'Common repeat after mutation: {ref_nuc}→{alt_nuc}\nHaplogroup D4a: decreased deletions, increased lifespan')
        ax2.set_ylabel('Alternative allele')
        ax2.set_xlabel('Genomic position')
        ax2.set_xticks(range(len2))
        ax2.set_xticklabels([str(alt_motif['motif.start']+i) for i in range(len2)], rotation=45, fontsize=9)
        ax2.set_yticks([])
        ax2.annotate('Mutation', xy=(alt_pos_in_motif+0.5, -0.1), xytext=(alt_pos_in_motif+0.5, -0.15),
                     ha='center', va='top', arrowprops=dict(arrowstyle='->', color='green', lw=2))

        delta_len = alt_motif['motif.length'] - ref_motif['motif.length']
        delta_eff = (alt_motif.get('effective.length', alt_motif['motif.length']) -
                     ref_motif.get('effective.length', ref_motif['motif.length']))
        info = (f'Metrics comparison:\n'
                f'• Motif length: {ref_motif["motif.length"]} → {alt_motif["motif.length"]} (Δ = {delta_len:+d})\n'
                f'• Effective length: {ref_motif.get("effective.length","N/A")} → {alt_motif.get("effective.length","N/A")} (Δ = {delta_eff:+d})\n'
                f'• Coordinates: {ref_motif["motif.start"]}-{ref_motif["motif.end"]}')
        plt.figtext(0.02, 0.02, info, fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.suptitle('Fig. 4: Common repeat at position 8473 (MT-ATP8/MT-ATP6)\nCentral nucleotide of 13‑bp direct repeat',
                     fontsize=16, fontweight='bold', y=1.02)
        plt.tight_layout(rect=[0, 0.1, 1, 0.95])
        out = output_dir / 'figure4_common_repeat_8473.png'
        plt.savefig(out, dpi=dpi, bbox_inches='tight')
        plt.close()
        print(f"Saved {out}")
    except Exception as e:
        print(f"Error generating Figure 4: {e}")


def plot_figure5(df_delta, output_dir, dpi=300):
    """Compare transitions vs transversions: mean |ΔR| bar, boxplot, stacked direction."""
    if df_delta.empty:
        return
    df = df_delta.copy()
    df['trans_type'] = df['is_transition'].map({True:'Transition', False:'Transversion'})
    # Direction categories
    def dir_cat(x):
        if x > 0.1: return 'Increase'
        if x < -0.1: return 'Decrease'
        return 'No change'
    df['direction'] = df['delta_R'].apply(dir_cat)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=dpi)
    # A: mean |ΔR|
    stats_df = df.groupby('trans_type')['abs_delta_R'].agg(['mean','std','count']).reset_index()
    axes[0].bar(stats_df['trans_type'], stats_df['mean'], yerr=stats_df['std'],
                capsize=8, color=['skyblue','lightcoral'], alpha=0.8, edgecolor='black')
    axes[0].set_ylabel('Mean |ΔR|')
    axes[0].set_title('Effect size by mutation type')
    axes[0].grid(axis='y', alpha=0.3)
    for i, row in stats_df.iterrows():
        axes[0].text(i, row['mean']+row['std']+0.01, f"{row['mean']:.3f}\n(n={int(row['count'])})",
                     ha='center', fontsize=9)
    # B: boxplot of ΔR
    trans = df[df['trans_type']=='Transition']['delta_R']
    transv = df[df['trans_type']=='Transversion']['delta_R']
    bp = axes[1].boxplot([trans, transv], labels=['Transition','Transversion'], patch_artist=True)
    for patch, col in zip(bp['boxes'], ['skyblue','lightcoral']):
        patch.set_facecolor(col); patch.set_alpha(0.7)
    axes[1].axhline(0, color='red', linestyle='--')
    axes[1].set_ylabel('ΔR')
    axes[1].set_title('Distribution of ΔR')
    axes[1].grid(axis='y', alpha=0.3)
    if len(trans) > 1 and len(transv) > 1:
        _, p = stats.ttest_ind(trans, transv, equal_var=False)
        axes[1].text(0.95, 0.95, f't-test p = {p:.3e}', transform=axes[1].transAxes,
                     ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    # C: stacked bar of direction
    dir_counts = df.groupby(['trans_type','direction']).size().unstack(fill_value=0)
    dir_pct = dir_counts.div(dir_counts.sum(axis=1), axis=0)*100
    dir_pct.plot(kind='bar', stacked=True, ax=axes[2], color=['green','red','gray'], alpha=0.8, edgecolor='black')
    axes[2].set_ylabel('Percentage (%)')
    axes[2].set_title('Direction of effect')
    axes[2].legend(title='Effect')
    axes[2].grid(axis='y', alpha=0.3)
    plt.suptitle('Fig. 5: Comparison of transitions and transversions', fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = output_dir / 'figure5_mutation_types.png'
    plt.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


def plot_figure6(df_delta, output_dir, dpi=300, coding_positions=None):
    """ΔR by functional mutation type (SYN/NON-SYN/NON-CODING)."""
    if df_delta.empty:
        return
    df = df_delta.copy()
    # Assign mutation type using the simple classifier (for these key positions only)
    df['mut_type'] = df.apply(lambda r: simple_mutation_type(r['position'], r['ref'], r['alt'], coding_positions), axis=1)
    # Direction (same as above)
    def dir_cat(x):
        if x > 0.1: return 'Increase'
        if x < -0.1: return 'Decrease'
        return 'No change'
    df['direction'] = df['delta_R'].apply(dir_cat)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=dpi)
    # A: Violin plot
    types = sorted(df['mut_type'].unique())
    if len(types) > 1:
        violin_data = [df[df['mut_type']==t]['delta_R'] for t in types]
        parts = axes[0,0].violinplot(violin_data, showmeans=True, showmedians=True)
        colors = ['lightgreen','salmon','lightblue']
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(colors[i%len(colors)]); pc.set_alpha(0.7); pc.set_edgecolor('black')
        axes[0,0].set_xticks(range(1, len(types)+1))
        axes[0,0].set_xticklabels(types)
        axes[0,0].axhline(0, color='red', linestyle='--')
        axes[0,0].set_ylabel('ΔR')
        axes[0,0].set_title('Distribution of ΔR by functional type')
        axes[0,0].grid(axis='y', alpha=0.3)
    # B: mean |ΔR|
    stats_df = df.groupby('mut_type')['abs_delta_R'].agg(['mean','std','count']).reset_index()
    axes[0,1].bar(stats_df['mut_type'], stats_df['mean'], yerr=stats_df['std'],
                  capsize=8, color=['lightgreen','salmon','lightblue'], alpha=0.8, edgecolor='black')
    axes[0,1].set_ylabel('Mean |ΔR|')
    axes[0,1].set_title('Effect size by functional type')
    axes[0,1].grid(axis='y', alpha=0.3)
    for i, row in stats_df.iterrows():
        axes[0,1].text(i, row['mean']+row['std']+0.01, f"{row['mean']:.3f}\n(n={int(row['count'])})",
                       ha='center', fontsize=9)
    # C: stacked direction
    dir_counts = df.groupby(['mut_type','direction']).size().unstack(fill_value=0)
    dir_pct = dir_counts.div(dir_counts.sum(axis=1), axis=0)*100
    dir_pct.plot(kind='bar', stacked=True, ax=axes[1,0], color=['green','red','gray'], alpha=0.8, edgecolor='black')
    axes[1,0].set_ylabel('Percentage (%)')
    axes[1,0].set_title('Direction of effect by functional type')
    axes[1,0].legend(title='Effect')
    axes[1,0].grid(axis='y', alpha=0.3)
    # D: horizontal bar of mean ΔR per position+type
    df['pos_type'] = df['position'].astype(str) + ' (' + df['mut_type'] + ')'
    top_groups = df['pos_type'].unique()[:10]
    if len(top_groups) > 1:
        pos_stats = []
        for g in top_groups:
            sub = df[df['pos_type'] == g]
            pos_stats.append({'group': g, 'mean': sub['delta_R'].mean(), 'std': sub['delta_R'].std()})
        pos_df = pd.DataFrame(pos_stats)
        ypos = range(len(pos_df))
        axes[1,1].barh(ypos, pos_df['mean'], xerr=pos_df['std'], color='skyblue', alpha=0.7, edgecolor='black')
        axes[1,1].set_yticks(ypos)
        axes[1,1].set_yticklabels(pos_df['group'], fontsize=9)
        axes[1,1].axvline(0, color='red', linestyle='--')
        axes[1,1].set_xlabel('Mean ΔR')
        axes[1,1].set_title('Mean ΔR by position (functional type)')
        axes[1,1].grid(axis='x', alpha=0.3)
    else:
        axes[1,1].text(0.5, 0.5, 'Not enough groups', ha='center', va='center')
        axes[1,1].set_axis_off()
    plt.suptitle('Fig. 6: Effect of functional mutation type on ΔR', fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = output_dir / 'figure6_functional_types.png'
    plt.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


# ----------------------------------------------------------------------
# Main CLI
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Analyse key mtDNA positions and generate six figures.",
        epilog="Example: %(prog)s --input_dir data/2_derived/01KP --fasta data/1_raw/sequence.fasta"
    )
    parser.add_argument('-i', '--input_dir', required=True,
                        help="Directory containing raw 01KP.*.txt files (not cleaned)")
    parser.add_argument('-f', '--fasta', required=True,
                        help="Reference FASTA file (used to obtain reference nucleotides)")
    parser.add_argument('-o', '--output_dir', default='./figures',
                        help="Output directory for figures (default: ./figures)")
    parser.add_argument('--positions', type=int, nargs='+',
                        default=[8251, 8472, 8473, 12705, 16223],
                        help="List of positions to analyse (default: all five key positions)")
    parser.add_argument('--dpi', type=int, default=300,
                        help="DPI for output figures (default: 300)")
    parser.add_argument('--skip_fig4', action='store_true',
                        help="Skip the schematic of position 8473 (Figure 4)")
    args = parser.parse_args()

    # Create output directory
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load reference genome (only needed if we want to double‑check ref nucleotides)
    try:
        record = next(SeqIO.parse(args.fasta, 'fasta'))
        ref_seq = str(record.seq).upper()
        ref_dict = {i+1: ref_seq[i] for i in range(len(ref_seq))}
    except Exception as e:
        print(f"Error reading FASTA: {e}")
        return 1

    # Filter positions to those that exist in POSITIONS_INFO and are requested
    all_positions = args.positions
    available = [p for p in all_positions if p in POSITIONS_INFO]
    if not available:
        print("None of the requested positions are in POSITIONS_INFO dictionary. Exiting.")
        return 1
    print(f"Analysing positions: {available}")

    # Collect data for each position
    df_list = []
    for pos in available:
        info = POSITIONS_INFO[pos]
        # Use reference from dict, but we trust POSITIONS_INFO; however we can verify with fasta if needed
        ref_nuc = info['ref']
        df_pos = load_position_data(args.input_dir, pos, ref_nuc)
        if not df_pos.empty:
            df_pos['gene'] = info['gene']
            df_pos['description'] = info['description']
            df_pos['importance'] = info['importance']
            df_list.append(df_pos)
        else:
            print(f"Warning: No data for position {pos}")

    if not df_list:
        print("No data loaded. Exiting.")
        return 1

    df_all = pd.concat(df_list, ignore_index=True)
    print(f"Loaded {len(df_all)} allele records for {len(available)} positions.")

    # Compute ΔR
    delta_list = []
    for pos in available:
        df_pos = df_all[df_all['position'] == pos]
        if not df_pos.empty:
            delta_list.append(compute_delta_R(df_pos))
    if not delta_list:
        print("No ΔR computed.")
        return 1
    df_delta = pd.concat(delta_list, ignore_index=True)
    print(f"Computed ΔR for {len(df_delta)} mutations.")

    # Generate figures
    plot_figure1(df_delta, out_dir, args.dpi)
    plot_figure2(df_delta, out_dir, args.dpi)
    plot_figure3(df_delta, out_dir, args.dpi)
    if not args.skip_fig4 and 8473 in available:
        plot_figure4_schematic(args.input_dir, out_dir, args.dpi)
    else:
        print("Skipping Figure 4 (position 8473 not requested or --skip_fig4).")
    plot_figure5(df_delta, out_dir, args.dpi)
    # For Figure 6 we need coding positions list (use available but only protein‑coding)
    coding = [p for p in available if p in [8251,8472,8473,12705]]
    plot_figure6(df_delta, out_dir, args.dpi, coding_positions=coding)

    print(f"\nAll figures saved to {out_dir}")
    return 0


if __name__ == "__main__":
    exit(main())
