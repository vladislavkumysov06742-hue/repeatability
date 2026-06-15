#!/usr/bin/env python3
"""
01C.KP.TryRepeatabilityMetrics.py

Calculate and compare repeatability for Ref and Alt alleles of specified mutations.
Generates statistics, plots, and saves results.

Author: Vladislav Gadzhiev
Date: 12/12/2025

Usage examples:
    python 01C_analyze.py --fasta ../data/1_raw/sequence.fasta --input_dir ../data/2_derived/01KP
    python 01C_analyze.py -f data.fasta -i repeats -o results --positions 8251 8473 12705
    python 01C_analyze.py --help
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import SeqIO, Seq
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


# ------------------------------
# 1. GENOME ANNOTATION (mitochondrial genes)
# ------------------------------
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
        return f"ERROR: ref mismatch (expected {seq[pos-1]}, got {ref})"
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
# 2. REPEATABILITY METRICS (using pre‑computed files)
# ------------------------------
def repeatability_perfect_longest(pos, nuc, input_dir, major_arc_start=5798, major_arc_end=16568):
    file_path = Path(input_dir) / f"01KP.{pos}.{nuc}.txt"
    try:
        df = pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        return np.nan
    mask = (~((df['motif.start'] == df['repeat.start']) & (df['motif.end'] == df['repeat.end']))) & \
           (df['motif.start'] >= major_arc_start) & (df['motif.end'] <= major_arc_end) & \
           (df['repeat.start'] >= major_arc_start) & (df['repeat.end'] <= major_arc_end) & \
           (df['repeat.hamming.distance'] == 0) & \
           (df['motif.start'] < df['pos']) & (df['pos'] < df['motif.end'])
    filtered = df[mask]
    if filtered.empty:
        return np.nan
    return filtered['motif.length'].max()


def repeatability_degraded_longest(pos, nuc, input_dir, max_degradation_pct=20,
                                   major_arc_start=5798, major_arc_end=16568):
    file_path = Path(input_dir) / f"01KP.{pos}.{nuc}.txt"
    try:
        df = pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        return np.nan
    df['degradation'] = ((df['motif.length'] - df['effective.length']) / df['motif.length']) * 100
    mask = (~((df['motif.start'] == df['repeat.start']) & (df['motif.end'] == df['repeat.end']))) & \
           (df['motif.start'] >= major_arc_start) & (df['motif.end'] <= major_arc_end) & \
           (df['repeat.start'] >= major_arc_start) & (df['repeat.end'] <= major_arc_end) & \
           (df['motif.start'] < df['pos']) & (df['pos'] < df['motif.end']) & \
           (df['degradation'] <= max_degradation_pct)
    filtered = df[mask]
    if filtered.empty:
        return np.nan
    return filtered['motif.length'].max()


# ------------------------------
# 3. MAIN ANALYSIS FUNCTION
# ------------------------------
def run_analysis(fasta_path, input_dir, output_dir, positions=None, max_degradation=20,
                 plot=True, save_plot=True):
    # Load reference
    record = next(SeqIO.parse(fasta_path, "fasta"))
    seq = str(record.seq).upper()
    print(f"Genome length: {len(seq)}")

    # Determine positions to analyze
    if positions is None:
        # Use all positions for which at least one repeat file exists
        files = list(Path(input_dir).glob("01KP.*.txt"))
        positions = sorted({int(f.stem.split('.')[1]) for f in files})
        print(f"Found {len(positions)} positions with repeat files.")
    else:
        positions = sorted(set(positions))
        print(f"Analyzing {len(positions)} user‑provided positions.")

    # Generate mutation table for these positions
    mutations = []
    for pos in positions:
        ref = seq[pos-1]
        if ref not in {'A','T','G','C'}:
            continue
        gene_info = get_gene_at_position(pos)
        gene_name = gene_info['name'] if gene_info else 'NON-CODING'
        gene_type = gene_info['type'] if gene_info else 'NON-CODING'
        for alt in ['A','T','G','C']:
            if alt == ref:
                continue
            if gene_type == 'protein_coding':
                mut_type = is_synonymous_mutation(seq, pos, ref, alt)
            else:
                mut_type = 'NON-CODING'
            mutations.append({
                'position': pos,
                'ref': ref,
                'alt': alt,
                'gene': gene_name,
                'mutation_type': mut_type,
            })
    mutations_df = pd.DataFrame(mutations)

    # Compute repeatability metrics
    results = []
    for _, row in mutations_df.iterrows():
        pos = row['position']
        ref = row['ref']
        alt = row['alt']
        ref_perf = repeatability_perfect_longest(pos, ref, input_dir)
        alt_perf = repeatability_perfect_longest(pos, alt, input_dir)
        ref_degr = repeatability_degraded_longest(pos, ref, input_dir, max_degradation)
        alt_degr = repeatability_degraded_longest(pos, alt, input_dir, max_degradation)

        results.append({
            'position': pos,
            'ref': ref,
            'alt': alt,
            'gene': row['gene'],
            'mutation_type': row['mutation_type'],
            'ref_perfect': ref_perf,
            'alt_perfect': alt_perf,
            'ref_degraded': ref_degr,
            'alt_degraded': alt_degr,
            'perfect_diff': alt_perf - ref_perf if not (np.isnan(alt_perf) or np.isnan(ref_perf)) else np.nan,
            'degraded_diff': alt_degr - ref_degr if not (np.isnan(alt_degr) or np.isnan(ref_degr)) else np.nan,
        })
    results_df = pd.DataFrame(results)

    # Create output directory
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Save data
    mutations_df.to_csv(out_path / 'mutations.csv', index=False)
    results_df.to_csv(out_path / 'delta_repeatability.csv', index=False)

    # Pattern analysis (optional)
    patterns = []
    for pos in positions:
        ref = seq[pos-1]
        start = max(0, pos-4)
        end = min(len(seq), pos+2)
        context = seq[start:end]
        gc = (context.count('G') + context.count('C')) / len(context) * 100
        for alt in ['A','T','G','C']:
            if alt == ref:
                continue
            mut_type = is_synonymous_mutation(seq, pos, ref, alt)
            trans = (ref, alt) in [('A','G'),('G','A'),('C','T'),('T','C')]
            patterns.append({
                'position': pos,
                'ref': ref,
                'alt': alt,
                'context': context,
                'gc_content': gc,
                'mutation_type': mut_type,
                'is_transition': trans,
            })
    patterns_df = pd.DataFrame(patterns)
    patterns_df.to_csv(out_path / 'mutation_patterns.csv', index=False)

    # Save Excel (optional)
    try:
        with pd.ExcelWriter(out_path / 'repeatability_analysis.xlsx') as writer:
            mutations_df.to_excel(writer, sheet_name='Mutations', index=False)
            results_df.to_excel(writer, sheet_name='DeltaR', index=False)
            patterns_df.to_excel(writer, sheet_name='Patterns', index=False)
        print(f"Excel file saved: {out_path / 'repeatability_analysis.xlsx'}")
    except Exception as e:
        print(f"Could not save Excel file: {e}")

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"Total mutations analyzed: {len(results_df)}")
    if 'perfect_diff' in results_df.columns:
        diff = results_df['perfect_diff'].dropna()
        if len(diff) > 0:
            print(f"ΔR (perfect) mean: {diff.mean():.3f}  median: {diff.median():.3f}")
            print(f"ΔR (perfect) std: {diff.std():.3f}")
            inc = (diff > 0).sum()
            dec = (diff < 0).sum()
            print(f"Increased: {inc} ({inc/len(diff)*100:.1f}%)  Decreased: {dec} ({dec/len(diff)*100:.1f}%)")
    # Mutation types
    print("\nBy mutation type:")
    for mt in results_df['mutation_type'].unique():
        sub = results_df[results_df['mutation_type'] == mt]
        print(f"  {mt}: n={len(sub)}  ΔR_mean={sub['perfect_diff'].mean():+.3f}")

    # Plotting
    if not plot:
        return

    # Create a multi‑panel figure (similar to notebook)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Repeatability changes by mutation type', fontsize=16)

    # Plot A: barplot of mean ΔR (perfect)
    mut_types = sorted(results_df['mutation_type'].unique())
    perfect_means = [results_df[results_df['mutation_type'] == t]['perfect_diff'].mean() for t in mut_types]
    bars1 = axes[0,0].bar(range(len(mut_types)), perfect_means, color='skyblue')
    axes[0,0].set_xticks(range(len(mut_types)))
    axes[0,0].set_xticklabels(mut_types, rotation=45)
    axes[0,0].set_ylabel('Mean ΔR (perfect)')
    axes[0,0].set_title('Ideal repeats')
    axes[0,0].axhline(0, color='r', linestyle='--')
    for i, v in enumerate(perfect_means):
        axes[0,0].text(i, v + 0.05*np.sign(v), f'{v:+.1f}', ha='center')

    # Plot B: barplot of mean ΔR (degraded)
    degraded_means = [results_df[results_df['mutation_type'] == t]['degraded_diff'].mean() for t in mut_types]
    bars2 = axes[0,1].bar(range(len(mut_types)), degraded_means, color='lightgreen')
    axes[0,1].set_xticks(range(len(mut_types)))
    axes[0,1].set_xticklabels(mut_types, rotation=45)
    axes[0,1].set_ylabel('Mean ΔR (degraded)')
    axes[0,1].set_title(f'Degraded repeats (≤{max_degradation}%)')
    axes[0,1].axhline(0, color='r', linestyle='--')
    for i, v in enumerate(degraded_means):
        axes[0,1].text(i, v + 0.05*np.sign(v), f'{v:+.1f}', ha='center')

    # Plot C: boxplot for SYN, NON-SYN, NON-CODING (perfect)
    syn = results_df[results_df['mutation_type'] == 'SYN']['perfect_diff'].dropna()
    nonsyn = results_df[results_df['mutation_type'] == 'NON-SYN']['perfect_diff'].dropna()
    noncoding = results_df[results_df['mutation_type'] == 'NON-CODING']['perfect_diff'].dropna()
    axes[1,0].boxplot([syn, nonsyn, noncoding], labels=['SYN','NON-SYN','NON-CODING'])
    axes[1,0].axhline(0, color='r', linestyle='--')
    axes[1,0].set_ylabel('ΔR (perfect)')
    axes[1,0].set_title('Distribution by functional class')
    # T‑test between SYN and NON‑SYN if enough data
    if len(syn) > 1 and len(nonsyn) > 1:
        tstat, pval = stats.ttest_ind(syn, nonsyn, nan_policy='omit')
        axes[1,0].text(0.5, 0.95, f'p = {pval:.3e}', transform=axes[1,0].transAxes,
                       ha='center', va='top', fontsize=9)

    # Plot D: pie chart of change direction (increase/decrease/no change)
    change_counts = results_df['perfect_change'].value_counts().to_dict()
    # Add missing categories
    for cat in ['INCREASE','DECREASE','NO_CHANGE']:
        if cat not in change_counts:
            change_counts[cat] = 0
    labels = list(change_counts.keys())
    sizes = list(change_counts.values())
    colors_pie = ['green','red','gray']
    axes[1,1].pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%')
    axes[1,1].set_title('Effect on perfect repeat length')

    plt.tight_layout()
    if save_plot:
        plot_path = out_path / 'repeatability_analysis_figure.png'
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved: {plot_path}")
    plt.show()


# ------------------------------
# 4. COMMAND LINE INTERFACE
# ------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Analyze repeatability changes for all possible nucleotide substitutions.",
        epilog="Example: %(prog)s --fasta ../data/1_raw/sequence.fasta --input_dir ../data/2_derived/01KP --output results"
    )
    parser.add_argument('-f', '--fasta', required=True, help="Reference FASTA file")
    parser.add_argument('-i', '--input_dir', required=True, help="Directory containing 01KP.*.txt repeat files")
    parser.add_argument('-o', '--output', default='repeatability_analysis', help="Output directory (default: repeatability_analysis)")
    parser.add_argument('-p', '--positions', type=int, nargs='+', help="List of positions (1‑based). If not given, use all positions with repeat files.")
    parser.add_argument('--max_degradation', type=int, default=20, help="Max degradation percent for degraded metric (default: 20)")
    parser.add_argument('--no_plot', action='store_true', help="Do not generate plots (only save data)")
    parser.add_argument('--save_plot', action='store_true', default=True, help="Save plot to file (default: True)")
    args = parser.parse_args()

    run_analysis(
        fasta_path=args.fasta,
        input_dir=args.input_dir,
        output_dir=args.output,
        positions=args.positions,
        max_degradation=args.max_degradation,
        plot=not args.no_plot,
        save_plot=args.save_plot
    )


if __name__ == "__main__":
    main()
