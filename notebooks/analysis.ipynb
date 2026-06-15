#!/usr/bin/env python3
"""
analysis.py

Comprehensive analysis of repeatability files: compute ΔR, generate statistical reports,
produce publication‑ready figures (including separate heatmaps per reference nucleotide,
and GC content analysis of mutation patterns).

Author: Vladislav Gadziev
Date: 15/12/2025

Usage examples:
    python analysis.py --input_dir data/2_derived/01KP/cleaned --fasta data/1_raw/sequence.fasta
    python analysis.py -i data/2_derived/01KP/cleaned -f mtDNA.fasta -o my_results --dpi 300
    python analysis.py --help
"""

import argparse
import sys
import re
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from Bio import SeqIO, Seq

warnings.filterwarnings('ignore')


# ------------------------------
# 1. Data loading and aggregation
# ------------------------------
def analyze_all_positions(data_dir, ref_seq_dict):
    """
    Read all 01KP.*.txt files, compute per‑allele summary statistics.
    Returns DataFrame and ref_nucleotides dictionary (from FASTA).
    """
    results = []
    data_path = Path(data_dir)
    file_pattern = re.compile(r"01KP\.(\d+)\.([ATGC])\.txt")
    files = list(data_path.glob("01KP.*.txt"))
    print(f"Found {len(files)} files.")

    ref_nucleotides = ref_seq_dict  # use provided reference

    for file_path in files:
        match = file_pattern.match(file_path.name)
        if not match:
            continue
        pos = int(match.group(1))
        nuc = match.group(2)

        try:
            df = pd.read_csv(file_path, sep='\t')
            if df.empty:
                continue

            # Basic stats
            mean_length = df['motif.length'].mean()
            max_length = df['motif.length'].max()
            median_length = df['motif.length'].median()

            if 'effective.length' in df.columns:
                mean_effective = df['effective.length'].mean()
                max_effective = df['effective.length'].max()
            else:
                mean_effective = mean_length
                max_effective = max_length

            # GC content of motifs
            gc_values = []
            for seq in df['motif.seq']:
                if isinstance(seq, str):
                    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
                    gc_values.append(gc)
            mean_gc = np.mean(gc_values) if gc_values else np.nan
            median_gc = np.median(gc_values) if gc_values else np.nan

            perfect_repeats = len(df[df['repeat.hamming.distance'] == 0]) if 'repeat.hamming.distance' in df.columns else np.nan

            total_bases = df['motif.length'].sum() if len(df) > 0 else 0
            avg_seq_len = total_bases / len(df) if len(df) > 0 else 0

            results.append({
                'position': pos,
                'nucleotide': nuc,
                'is_reference': (ref_nucleotides.get(pos) == nuc),
                'mean_length': mean_length,
                'max_length': max_length,
                'median_length': median_length,
                'mean_effective': mean_effective,
                'max_effective': max_effective,
                'mean_gc': mean_gc,
                'median_gc': median_gc,
                'perfect_repeats': perfect_repeats,
                'total_repeats': len(df),
                'percent_perfect': (perfect_repeats / len(df) * 100) if len(df) > 0 else 0,
                'total_bases': total_bases,
                'avg_sequence_length': avg_seq_len,
                'file_path': str(file_path)
            })
        except Exception as e:
            print(f"Error reading {file_path.name}: {e}")

    return pd.DataFrame(results), ref_nucleotides


# ------------------------------
# 2. ΔR calculation
# ------------------------------
def is_transition(ref, alt):
    return (ref, alt) in {('A','G'), ('G','A'), ('C','T'), ('T','C')}


def calculate_delta_R(df_all, ref_nucleotides):
    """Compute ΔR for every alternative allele relative to the reference."""
    delta_data = []
    for pos in df_all['position'].unique():
        pos_data = df_all[df_all['position'] == pos]
        ref_nuc = ref_nucleotides.get(pos)
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


# ------------------------------
# 3. GC content analysis of mutation contexts
# ------------------------------
def get_mutation_context(seq, pos, left=4, right=2):
    """Extract context around position (left bases upstream, right bases downstream).
    Default left=4, right=2 as in original notebook."""
    start = max(0, pos - left - 1)   # -1 for 0‑based? In notebook: start = max(0, pos-4)
    # Actually in notebook: start = max(0, pos-4) (pos is 1‑based, so pos-4 gives 0‑based start)
    # Let's replicate exactly: start = max(0, pos-4)   (since pos 1‑based)
    # End = min(len(seq), pos+2)   (pos+2 is exclusive in slice)
    start = max(0, pos - left)
    end = min(len(seq), pos + right)
    context = seq[start:end]
    return context


def analyze_mutation_patterns(seq, df_delta, positions=None):
    """
    Analyze mutation patterns: GC content of context, relationship with mutation type.
    Similar to the 'analyze_mutation_patterns' function in the original notebook.
    """
    if positions is None:
        positions = df_delta['position'].unique()
    patterns = []
    for pos in positions:
        ref = seq[pos-1]
        start = max(0, pos-4)
        end = min(len(seq), pos+2)
        context = seq[start:end]
        gc_content = (context.count('G') + context.count('C')) / len(context) * 100 if len(context)>0 else np.nan
        for alt in ['A','T','G','C']:
            if alt == ref:
                continue
            # Determine mutation type (synonymous/non-synonymous etc.) using imported functions
            # We need a simple helper; for simplicity we can reuse is_synonymous_mutation logic.
            # Here I'll use the original function from the notebook (copied below)
            mut_type = is_synonymous_mutation(seq, pos, ref, alt)
            patterns.append({
                'position': pos,
                'ref': ref,
                'alt': alt,
                'context': context,
                'gc_content': gc_content,
                'mutation_type': mut_type,
                'is_transition': is_transition(ref, alt),
                'is_transversion': not is_transition(ref, alt)
            })
    return pd.DataFrame(patterns)


# Synonymous mutation detection (from original notebook)
mitochondrial_genes = [
    {'name': 'MT-ND1',  'start': 3307,  'end': 4262,  'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-ND2',  'start': 4470,  'end': 5511,  'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-CO1',  'start': 5904,  'end': 7445,  'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-CO2',  'start': 7586,  'end': 8269,  'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-ATP8', 'start': 8366,  'end': 8572,  'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-ATP6', 'start': 8527,  'end': 9207,  'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-CO3',  'start': 9207,  'end': 9990,  'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-ND3',  'start': 10059, 'end': 10404, 'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-ND4L', 'start': 10470, 'end': 10766, 'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-ND4',  'start': 10760, 'end': 12137, 'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-ND5',  'start': 12337, 'end': 14148, 'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-ND6',  'start': 14149, 'end': 14673, 'type': 'protein_coding', 'strand': '-'},
    {'name': 'MT-CYB',  'start': 14747, 'end': 15887, 'type': 'protein_coding', 'strand': '+'},
    {'name': 'MT-RNR1', 'start': 648,   'end': 1601,  'type': 'rRNA', 'strand': '+'},
    {'name': 'MT-RNR2', 'start': 1671,  'end': 3229,  'type': 'rRNA', 'strand': '+'},
    {'name': 'MT-TF',   'start': 577,   'end': 647,   'type': 'tRNA', 'strand': '+'},
    {'name': 'MT-TV',   'start': 1602,  'end': 1670,  'type': 'tRNA', 'strand': '+'},
]

mitochondrial_genetic_code = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': '*', 'AGG': '*',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'W', 'TGG': 'W',
}
stop_codons = {'TAA', 'TAG', 'AGA', 'AGG'}

def get_gene_at_position(pos):
    for gene in mitochondrial_genes:
        if gene['start'] <= pos <= gene['end']:
            return gene
    return None

def get_codon_at_position(seq, pos, gene_info):
    if gene_info['type'] != 'protein_coding':
        return None, None, None
    if gene_info['strand'] == '+':
        gene_pos = pos - gene_info['start'] + 1
    else:
        gene_pos = gene_info['end'] - pos + 1
    codon_num = (gene_pos + 2) // 3
    if gene_info['strand'] == '+':
        codon_start = gene_info['start'] + (codon_num - 1) * 3
    else:
        codon_start = gene_info['end'] - (codon_num - 1) * 3 - 2
    if gene_info['strand'] == '+':
        codon = seq[codon_start-1:codon_start+2]
    else:
        codon_seq = seq[codon_start-1:codon_start+2]
        codon = str(Seq.Seq(codon_seq).reverse_complement())
    codon_pos = (gene_pos - 1) % 3
    return codon, codon_pos, codon_num

def codon_to_amino_acid(codon):
    return mitochondrial_genetic_code.get(codon.upper(), 'X')

def is_synonymous_mutation(seq, pos, ref, alt):
    if seq[pos-1] != ref:
        return f"ERROR: ref mismatch"
    gene_info = get_gene_at_position(pos)
    if not gene_info or gene_info['type'] != 'protein_coding':
        return 'NON-CODING'
    original_codon, codon_pos, _ = get_codon_at_position(seq, pos, gene_info)
    if original_codon is None:
        return 'UNKNOWN'
    mutated_codon = list(original_codon)
    mutated_codon[codon_pos] = alt
    mutated_codon = ''.join(mutated_codon)
    original_aa = codon_to_amino_acid(original_codon)
    mutated_aa = codon_to_amino_acid(mutated_codon)
    if original_aa == 'X' or mutated_aa == 'X':
        return 'UNKNOWN'
    orig_stop = original_codon in stop_codons
    mut_stop = mutated_codon in stop_codons
    if orig_stop and not mut_stop:
        return 'STOP-LOSS'
    if not orig_stop and mut_stop:
        return 'STOP-GAIN'
    if original_aa == mutated_aa:
        return 'SYN'
    return 'NON-SYN'


# ------------------------------
# 4. Text report (extended)
# ------------------------------
def create_extended_text_report(df_delta, df_all, ref_nucleotides, patterns_df, output_path):
    """Write a comprehensive text report for machine learning and interpretation,
    including GC content analysis results."""
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("=" * 100 + "\n")
        f.write("FULL ΔR ANALYSIS REPORT FOR MACHINE LEARNING\n")
        f.write("=" * 100 + "\n\n")

        # 1. Dataset overview
        f.write("1. DATASET OVERVIEW\n")
        f.write("-" * 60 + "\n")
        f.write(f"Total positions in genome: {df_all['position'].nunique()}\n")
        f.write(f"Total alleles (position + nucleotide): {len(df_all)}\n")
        f.write(f"Unique nucleotides: {len(df_all['nucleotide'].unique())}\n")
        f.write(f"Average alleles per position: {len(df_all) / df_all['position'].nunique():.2f}\n")
        f.write(f"\nReference nucleotides distribution:\n")
        ref_counts = Counter(ref_nucleotides.values())
        for nuc in ['A','T','G','C']:
            cnt = ref_counts.get(nuc, 0)
            f.write(f"  {nuc}: {cnt} positions ({cnt/len(ref_nucleotides)*100:.1f}%)\n")

        # 2. ΔR statistics
        f.write("\n2. ΔR STATISTICS (CHANGE IN REPEAT LENGTH)\n")
        f.write("-" * 60 + "\n")
        if len(df_delta) > 0:
            stats_dict = {
                'mean': df_delta['delta_R'].mean(),
                'median': df_delta['delta_R'].median(),
                'std': df_delta['delta_R'].std(),
                'min': df_delta['delta_R'].min(),
                'max': df_delta['delta_R'].max(),
                'q25': df_delta['delta_R'].quantile(0.25),
                'q75': df_delta['delta_R'].quantile(0.75),
                'iqr': df_delta['delta_R'].quantile(0.75) - df_delta['delta_R'].quantile(0.25),
                'abs_mean': df_delta['abs_delta_R'].mean(),
                'abs_median': df_delta['abs_delta_R'].median()
            }
            for k, v in stats_dict.items():
                f.write(f"{k}: {v:.4f}\n")
            inc = len(df_delta[df_delta['delta_R'] > 0.1])
            dec = len(df_delta[df_delta['delta_R'] < -0.1])
            neut = len(df_delta) - inc - dec
            f.write(f"\nDirection of change (threshold ±0.1):\n")
            f.write(f"  Increase (>0.1): {inc} ({inc/len(df_delta)*100:.1f}%)\n")
            f.write(f"  Decrease (<-0.1): {dec} ({dec/len(df_delta)*100:.1f}%)\n")
            f.write(f"  Neutral: {neut} ({neut/len(df_delta)*100:.1f}%)\n")

        # 3. Mutation type analysis
        f.write("\n3. MUTATION TYPE ANALYSIS\n")
        f.write("-" * 60 + "\n")
        if len(df_delta) > 0:
            trans = df_delta[df_delta['is_transition']]
            transv = df_delta[~df_delta['is_transition']]
            f.write(f"Transitions: {len(trans)} ({len(trans)/len(df_delta)*100:.1f}%)\n")
            f.write(f"Transversions: {len(transv)} ({len(transv)/len(df_delta)*100:.1f}%)\n")
            if len(trans) > 0:
                f.write(f"\nTransitions:\n  mean ΔR = {trans['delta_R'].mean():.4f}, mean |ΔR| = {trans['abs_delta_R'].mean():.4f}, std = {trans['delta_R'].std():.4f}\n")
            if len(transv) > 0:
                f.write(f"\nTransversions:\n  mean ΔR = {transv['delta_R'].mean():.4f}, mean |ΔR| = {transv['abs_delta_R'].mean():.4f}, std = {transv['delta_R'].std():.4f}\n")

        # 4. Top positions
        f.write("\n4. TOP POSITIONS BY EFFECT SIZE\n")
        f.write("-" * 60 + "\n")
        if len(df_delta) > 0:
            pos_stats = df_delta.groupby('position').agg(
                mean_abs_delta=('abs_delta_R','mean'),
                mean_delta=('delta_R','mean'),
                count=('delta_R','count')
            ).round(4)
            top10 = pos_stats.nlargest(10, 'mean_abs_delta')
            for pos, row in top10.iterrows():
                f.write(f"  Position {pos}: |ΔR| = {row['mean_abs_delta']:.3f}, ΔR = {row['mean_delta']:.3f}, n = {int(row['count'])}\n")

        # 5. GC content analysis of mutation contexts
        f.write("\n5. GC CONTENT ANALYSIS OF MUTATION CONTEXTS\n")
        f.write("-" * 60 + "\n")
        if patterns_df is not None and not patterns_df.empty:
            f.write(f"Total context records: {len(patterns_df)}\n")
            f.write(f"Mean GC content: {patterns_df['gc_content'].mean():.2f}%\n")
            # Relation with mutation type
            syn = patterns_df[patterns_df['mutation_type'] == 'SYN']
            nonsyn = patterns_df[patterns_df['mutation_type'] == 'NON-SYN']
            if len(syn) > 0 and len(nonsyn) > 0:
                f.write(f"\nSYN mutations: mean GC context = {syn['gc_content'].mean():.2f}% (n={len(syn)})\n")
                f.write(f"NON-SYN mutations: mean GC context = {nonsyn['gc_content'].mean():.2f}% (n={len(nonsyn)})\n")
                # test
                if len(syn) > 5 and len(nonsyn) > 5:
                    _, p = stats.ttest_ind(syn['gc_content'], nonsyn['gc_content'], nan_policy='omit')
                    f.write(f"t-test p-value: {p:.4e}\n")
            # GC threshold analysis (as in original notebook)
            for thr in [30, 50, 70]:
                high_gc = patterns_df[patterns_df['gc_content'] >= thr]
                low_gc = patterns_df[patterns_df['gc_content'] < thr]
                if len(high_gc) > 0 and len(low_gc) > 0:
                    high_syn_pct = len(high_gc[high_gc['mutation_type'] == 'SYN']) / len(high_gc) * 100
                    low_syn_pct = len(low_gc[low_gc['mutation_type'] == 'SYN']) / len(low_gc) * 100
                    f.write(f"\nGC ≥ {thr}%: {high_syn_pct:.1f}% SYN (n={len(high_gc)})\n")
                    f.write(f"GC < {thr}%: {low_syn_pct:.1f}% SYN (n={len(low_gc)})\n")

        # 6. Correlation analysis
        f.write("\n6. CORRELATION ANALYSIS\n")
        f.write("-" * 60 + "\n")
        if len(df_delta) > 10:
            corr = df_delta[['delta_R','abs_delta_R','ref_length','delta_GC','delta_perfect']].corr()
            f.write("Pearson correlation matrix (upper triangle):\n")
            for i, col1 in enumerate(corr.columns):
                for j, col2 in enumerate(corr.columns):
                    if i < j:
                        f.write(f"  {col1} vs {col2}: {corr.iloc[i,j]:.4f}\n")
            # statistical significance for selected pairs
            f.write("\nStatistical significance (p-value):\n")
            pairs = [('ref_length','abs_delta_R'), ('delta_GC','delta_R'), ('ref_length','delta_R'), ('delta_perfect','delta_R')]
            for c1, c2 in pairs:
                valid = df_delta[[c1,c2]].dropna()
                if len(valid) > 2:
                    r, p = stats.pearsonr(valid[c1], valid[c2])
                    f.write(f"  {c1} vs {c2}: r={r:.4f}, p={p:.4e}\n")

        # 7. Outliers
        f.write("\n7. OUTLIERS (1.5×IQR rule)\n")
        f.write("-" * 60 + "\n")
        if len(df_delta) > 0:
            Q1 = df_delta['delta_R'].quantile(0.25)
            Q3 = df_delta['delta_R'].quantile(0.75)
            IQR = Q3 - Q1
            lower = Q1 - 1.5*IQR
            upper = Q3 + 1.5*IQR
            out = df_delta[(df_delta['delta_R'] < lower) | (df_delta['delta_R'] > upper)]
            f.write(f"Outliers: {len(out)} ({len(out)/len(df_delta)*100:.1f}%)\n")
            if len(out) > 0:
                f.write("Top 5 positive outliers:\n")
                for _, row in out.nlargest(5, 'delta_R').iterrows():
                    f.write(f"  pos {row['position']}, {row['mutation']}: ΔR={row['delta_R']:.3f}\n")
                f.write("Top 5 negative outliers:\n")
                for _, row in out.nsmallest(5, 'delta_R').iterrows():
                    f.write(f"  pos {row['position']}, {row['mutation']}: ΔR={row['delta_R']:.3f}\n")

        # 8. Recommendations
        f.write("\n8. RECOMMENDATIONS FOR ML MODELS\n")
        f.write("-" * 60 + "\n")
        f.write("Potential features:\n")
        f.write("  1. genomic position\n  2. reference nucleotide (one-hot)\n  3. alternative nucleotide\n")
        f.write("  4. transition/transversion flag\n  5. reference repeat length\n  6. reference GC content\n")
        f.write("  7. percent perfect repeats in reference\n  8. position normalised by genome length\n")
        f.write("  9. GC content of local context (from patterns analysis)\n")
        if df_delta.isnull().sum().sum() > 0:
            f.write(f"\nMissing values: {df_delta.isnull().sum().sum()}\n")
        f.write("\nPreprocessing steps:\n  1. Normalise numeric features\n  2. One‑hot encode categorical variables\n")
        f.write("  3. Handle outliers\n  4. Check for multicollinearity\n")
        f.write("\n" + "="*100 + "\nEND OF REPORT\n" + "="*100 + "\n")


# ------------------------------
# 5. Plotting functions
# ------------------------------
def plot_figure1(df_delta, output_dir, dpi):
    """Scatter plot (ref_length vs |ΔR|) + histogram of ΔR."""
    if len(df_delta) < 10:
        print("Not enough data for Figure 1.")
        return
    top_pos = df_delta.groupby('position').size().nlargest(100).index
    df_top = df_delta[df_delta['position'].isin(top_pos)]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=dpi)
    # scatter
    sc = axes[0].scatter(df_top['ref_length'], df_top['abs_delta_R'],
                         c=df_top['delta_GC'], cmap='coolwarm',
                         alpha=0.7, s=50, edgecolor='black', linewidth=0.5)
    axes[0].set_xlabel('Mean reference repeat length (bp)')
    axes[0].set_ylabel('|ΔR|')
    axes[0].set_title('|ΔR| vs reference repeat length (top 100 positions)')
    axes[0].grid(True, alpha=0.3)
    if df_top['ref_length'].nunique() > 1 and df_top['abs_delta_R'].nunique() > 1:
        slope, intercept, r_val, p_val, _ = stats.linregress(df_top['ref_length'], df_top['abs_delta_R'])
        x_line = np.array([df_top['ref_length'].min(), df_top['ref_length'].max()])
        axes[0].plot(x_line, intercept + slope*x_line, 'r--', label=f'R²={r_val**2:.3f}, p={p_val:.2e}')
        axes[0].legend()
    plt.colorbar(sc, ax=axes[0], label='ΔGC (%)')

    # histogram
    sns.histplot(df_delta['delta_R'], bins=50, kde=True, color='steelblue', edgecolor='black', ax=axes[1])
    axes[1].axvline(0, color='red', linestyle='--')
    axes[1].set_xlabel('ΔR')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Distribution of ΔR (all mutations)')
    axes[1].grid(True, alpha=0.3)
    mean_d = df_delta['delta_R'].mean()
    std_d = df_delta['delta_R'].std()
    axes[1].text(0.95, 0.95, f'μ = {mean_d:.3f}\nσ = {std_d:.3f}\nn = {len(df_delta)}',
                 transform=axes[1].transAxes, ha='right', va='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    plt.tight_layout()
    plt.savefig(output_dir / 'figure1_all_positions.png', dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Figure 1 saved: {output_dir / 'figure1_all_positions.png'}")


def plot_individual_heatmaps(df_delta, output_dir, dpi):
    """Generate separate heatmaps for each reference nucleotide (G, C, A, T)
    exactly as in the original notebook."""
    ref_list = ['G','C','A','T']
    possible_alts = {'G':['A','C','T'], 'C':['A','G','T'], 'A':['C','G','T'], 'T':['A','C','G']}
    for ref_nuc in ref_list:
        df_ref = df_delta[df_delta['ref'] == ref_nuc]
        if df_ref.empty:
            print(f"No data for reference {ref_nuc}, skipping heatmap.")
            continue
        alts = possible_alts[ref_nuc]
        # Top 20 positions by mean |ΔR| for this ref
        top_pos = df_ref.groupby('position')['abs_delta_R'].mean().nlargest(20).index
        df_heat = df_ref[df_ref['position'].isin(top_pos)]
        if df_heat.empty:
            print(f"Insufficient data for {ref_nuc} heatmap.")
            continue
        pivot = df_heat.pivot_table(index='position', columns='mutation', values='abs_delta_R', aggfunc='mean', fill_value=0)
        for alt in alts:
            col = f"{ref_nuc}→{alt}"
            if col not in pivot.columns:
                pivot[col] = 0
        pivot = pivot[[f"{ref_nuc}→{alt}" for alt in alts]]
        # Sort positions by mean |ΔR| descending
        pivot['mean'] = pivot.mean(axis=1)
        pivot = pivot.sort_values('mean', ascending=False).drop(columns='mean')

        fig, ax = plt.subplots(figsize=(10, 8), dpi=dpi)
        im = ax.imshow(pivot.values, cmap='YlOrRd', aspect='auto', vmin=0, vmax=max(1, pivot.values.max()))
        ax.set_xticks(range(pivot.shape[1]))
        ax.set_yticks(range(pivot.shape[0]))
        ax.set_xticklabels(pivot.columns, rotation=45, ha='right')
        ax.set_yticklabels(pivot.index)
        ax.set_xlabel('Mutation type')
        ax.set_ylabel('Position')
        ax.set_title(f'Heatmap: mean |ΔR| for reference {ref_nuc}')
        plt.colorbar(im, ax=ax, shrink=0.8, label='Mean |ΔR|')
        # Add value labels inside cells
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                val = pivot.iloc[i, j]
                if val > 0:
                    textcol = 'white' if val > pivot.values.max()/2 else 'black'
                    ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=9, color=textcol, weight='bold')
        plt.tight_layout()
        out_file = output_dir / f'figure3_heatmap_ref_{ref_nuc}.png'
        plt.savefig(out_file, dpi=dpi, bbox_inches='tight')
        plt.close()
        print(f"Heatmap saved: {out_file}")


def plot_figure5(df_delta, output_dir, dpi):
    """Mutation type analysis: barplot, boxplot, stacked bar, top10 counts."""
    if df_delta.empty:
        return
    df = df_delta.copy()
    df['transition_type'] = df['is_transition'].map({True:'Transition', False:'Transversion'})
    # categorize change direction
    def cat_delta(x):
        if x > 1.0: return 'Strong increase'
        if x > 0.1: return 'Moderate increase'
        if x > 0.01: return 'Weak increase'
        if abs(x) <= 0.01: return 'No change'
        if x > -0.1: return 'Weak decrease'
        if x > -1.0: return 'Moderate decrease'
        return 'Strong decrease'
    df['change_dir'] = df['delta_R'].apply(cat_delta)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=dpi)

    # A: mean |ΔR| by transition type
    type_stats = df.groupby('transition_type')['abs_delta_R'].agg(['mean','std','count']).reset_index()
    axes[0,0].bar(type_stats['transition_type'], type_stats['mean'], yerr=type_stats['std'],
                  capsize=8, color=['skyblue','lightcoral'], alpha=0.8, edgecolor='black')
    axes[0,0].set_ylabel('Mean |ΔR|')
    axes[0,0].set_title('Effect size by mutation type')
    axes[0,0].grid(axis='y', alpha=0.3)
    for i, row in type_stats.iterrows():
        axes[0,0].text(i, row['mean']+row['std']+0.01, f"{row['mean']:.3f}\n(n={int(row['count'])})",
                       ha='center', fontsize=9)

    # B: boxplot of ΔR
    trans_data = df[df['transition_type']=='Transition']['delta_R']
    transv_data = df[df['transition_type']=='Transversion']['delta_R']
    bp = axes[0,1].boxplot([trans_data, transv_data], labels=['Transition','Transversion'], patch_artist=True)
    for patch, col in zip(bp['boxes'], ['skyblue','lightcoral']):
        patch.set_facecolor(col)
        patch.set_alpha(0.7)
    axes[0,1].axhline(0, color='red', linestyle='--')
    axes[0,1].set_ylabel('ΔR')
    axes[0,1].set_title('Distribution of ΔR')
    axes[0,1].grid(axis='y', alpha=0.3)

    # C: stacked bar (direction by type)
    dir_counts = df.groupby(['transition_type','change_dir']).size().unstack(fill_value=0)
    order = ['Strong increase','Moderate increase','Weak increase','No change','Weak decrease','Moderate decrease','Strong decrease']
    dir_counts = dir_counts[[c for c in order if c in dir_counts.columns]]
    dir_pct = dir_counts.div(dir_counts.sum(axis=1), axis=0) * 100
    colors_map = {'Strong increase':'#1b5e20','Moderate increase':'#4caf50','Weak increase':'#81c784',
                  'No change':'#bdbdbd','Weak decrease':'#ffcc80','Moderate decrease':'#ff9800','Strong decrease':'#e65100'}
    colors_plot = [colors_map[c] for c in dir_pct.columns]
    bottom = np.zeros(len(dir_pct))
    for i, col in enumerate(dir_pct.columns):
        axes[1,0].bar(dir_pct.index, dir_pct[col], bottom=bottom, color=colors_plot[i], alpha=0.8, edgecolor='black', width=0.6)
        bottom += dir_pct[col].values
    axes[1,0].set_ylabel('Percentage (%)')
    axes[1,0].set_title('Direction of effect by mutation type')
    axes[1,0].set_ylim(0, 105)
    axes[1,0].grid(axis='y', alpha=0.3)
    # add percentage labels inside bars
    for i, (idx, row) in enumerate(dir_pct.iterrows()):
        cum = 0
        for j, col in enumerate(dir_pct.columns):
            val = row[col]
            if val > 5:
                ypos = cum + val/2
                axes[1,0].text(i, ypos, f'{val:.1f}%', ha='center', va='center', fontsize=8,
                               color='white' if val > 20 else 'black')
            cum += val
    from matplotlib.patches import Rectangle
    leg = [Rectangle((0,0),1,1, facecolor=colors_map[c], edgecolor='black', alpha=0.8) for c in dir_pct.columns]
    axes[1,0].legend(leg, dir_pct.columns, title='Effect direction', loc='upper left', bbox_to_anchor=(1.02,1))

    # D: top10 mutations by frequency
    mut_counts = df['mutation'].value_counts().head(10)
    axes[1,1].barh(range(len(mut_counts)), mut_counts.values, color='mediumpurple', alpha=0.7, edgecolor='black')
    axes[1,1].set_yticks(range(len(mut_counts)))
    axes[1,1].set_yticklabels(mut_counts.index)
    axes[1,1].set_xlabel('Count')
    axes[1,1].set_title('Top 10 most frequent mutations')
    for i, v in enumerate(mut_counts.values):
        axes[1,1].text(v + max(mut_counts.values)*0.01, i, f'{v}', va='center', fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_dir / 'figure5_mutation_analysis.png', dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Figure 5 saved: {output_dir / 'figure5_mutation_analysis.png'}")


def plot_gc_content_analysis(patterns_df, output_dir, dpi):
    """Plot GC content distribution and stacked bar of mutation types by GC bins."""
    if patterns_df is None or patterns_df.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=dpi)
    # Histogram of GC content
    axes[0].hist(patterns_df['gc_content'].dropna(), bins=30, color='darkgreen', alpha=0.7, edgecolor='black')
    axes[0].set_xlabel('GC content of context (%)')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Distribution of GC content around mutated positions')
    axes[0].grid(alpha=0.3)
    # Mutation type by GC bins
    bins = np.arange(0, 101, 10)
    patterns_df['gc_bin'] = pd.cut(patterns_df['gc_content'], bins)
    type_by_bin = patterns_df.groupby(['gc_bin', 'mutation_type']).size().unstack(fill_value=0)
    type_by_bin_pct = type_by_bin.div(type_by_bin.sum(axis=1), axis=0) * 100
    type_by_bin_pct.plot(kind='bar', stacked=True, ax=axes[1], colormap='viridis', alpha=0.8, edgecolor='black')
    axes[1].set_xlabel('GC content bin (%)')
    axes[1].set_ylabel('Percentage')
    axes[1].set_title('Mutation type composition by GC content')
    axes[1].legend(title='Mutation type', bbox_to_anchor=(1.05, 1), loc='upper left')
    axes[1].tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.savefig(output_dir / 'gc_content_analysis.png', dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"GC content analysis plot saved: {output_dir / 'gc_content_analysis.png'}")


# ------------------------------
# 6. Main CLI
# ------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive analysis of repeatability: compute ΔR, generate statistics and figures.",
        epilog="Example: %(prog)s --input_dir data/2_derived/01KP/cleaned --fasta data/1_raw/sequence.fasta --outdir results"
    )
    parser.add_argument('-i', '--input_dir', required=True, help="Directory containing cleaned 01KP.*.txt files")
    parser.add_argument('-f', '--fasta', required=True, help="Reference FASTA file (used for reference nucleotides)")
    parser.add_argument('-o', '--outdir', default='./repeatability_analysis', help="Output directory (default: ./repeatability_analysis)")
    parser.add_argument('--dpi', type=int, default=300, help="DPI for figures (default: 300)")
    parser.add_argument('--skip_plots', action='store_true', help="Skip generating plots (only text report and CSV)")
    parser.add_argument('--skip_report', action='store_true', help="Skip text report (only CSV and plots)")
    args = parser.parse_args()

    # Create output directories
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    text_dir = out_dir / 'text_reports'
    if not args.skip_report:
        text_dir.mkdir(exist_ok=True)

    # Load reference genome
    try:
        record = next(SeqIO.parse(args.fasta, 'fasta'))
        ref_seq = str(record.seq).upper()
        ref_nuc_dict = {i+1: ref_seq[i] for i in range(len(ref_seq))}
        print(f"Reference genome length: {len(ref_seq)}")
    except Exception as e:
        print(f"Error reading FASTA: {e}")
        return 1

    # Load and aggregate all data
    print("Loading repeatability files...")
    df_all, ref_nucs = analyze_all_positions(args.input_dir, ref_nuc_dict)
    if df_all.empty:
        print("No data loaded. Exiting.")
        return 1
    print(f"Loaded {len(df_all)} allele records from {df_all['position'].nunique()} positions.")

    # Calculate ΔR
    df_delta = calculate_delta_R(df_all, ref_nucs)
    if df_delta.empty:
        print("No ΔR values computed (no alternative alleles?). Exiting.")
        return 1
    print(f"Computed ΔR for {len(df_delta)} mutations.")

    # Save raw data
    df_all.to_csv(out_dir / 'all_positions_data.csv', index=False)
    df_delta.to_csv(out_dir / 'delta_R_data.csv', index=False)
    print(f"Saved CSV files to {out_dir}")

    # GC content analysis of mutation contexts
    print("Analyzing mutation contexts (GC content)...")
    patterns_df = analyze_mutation_patterns(ref_seq, df_delta)
    if not patterns_df.empty:
        patterns_df.to_csv(out_dir / 'mutation_patterns.csv', index=False)
        print(f"Saved mutation patterns to {out_dir / 'mutation_patterns.csv'}")
    else:
        patterns_df = None

    # Text report (includes GC analysis)
    if not args.skip_report:
        report_path = text_dir / 'full_analysis_report.txt'
        create_extended_text_report(df_delta, df_all, ref_nucs, patterns_df, report_path)
        print(f"Text report saved: {report_path}")

    # Figures
    if not args.skip_plots:
        print("Generating figures...")
        plot_figure1(df_delta, out_dir, args.dpi)
        plot_individual_heatmaps(df_delta, out_dir, args.dpi)   # separate heatmaps per ref
        plot_figure5(df_delta, out_dir, args.dpi)
        if patterns_df is not None and not patterns_df.empty:
            plot_gc_content_analysis(patterns_df, out_dir, args.dpi)

    # Final console summary
    print("\n" + "="*60)
    print("FINAL SUMMARY")
    print("="*60)
    print(f"Total positions: {df_all['position'].nunique()}")
    print(f"Total mutations analyzed: {len(df_delta)}")
    print(f"Mean ΔR: {df_delta['delta_R'].mean():.4f}  median: {df_delta['delta_R'].median():.4f}")
    inc = (df_delta['delta_R'] > 0).sum()
    dec = (df_delta['delta_R'] < 0).sum()
    print(f"Increased: {inc} ({inc/len(df_delta)*100:.1f}%)  Decreased: {dec} ({dec/len(df_delta)*100:.1f}%)")
    print(f"Transitions: {(df_delta['is_transition']).sum()} ({(df_delta['is_transition']).sum()/len(df_delta)*100:.1f}%)")
    if patterns_df is not None and not patterns_df.empty:
        mean_gc = patterns_df['gc_content'].mean()
        print(f"Mean GC content of mutation contexts: {mean_gc:.2f}%")
    print(f"All results saved in {out_dir}")
    return 0


if __name__ == "__main__":
    exit(main())
