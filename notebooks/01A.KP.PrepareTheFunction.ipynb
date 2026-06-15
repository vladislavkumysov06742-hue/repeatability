#!/usr/bin/env python3
"""
01A.KP.PrepareTheFunction.py

Calculate greedy repeatability for a given allele in a given position.

Author: Vladislav Gadzhiev
Date: 10/12/2025

Usage examples:
    python 01A_prepare_function.py --fasta ../data/1_raw/sequence.fasta --pos 8473 --nuc T
    python 01A_prepare_function.py -f ../data/1_raw/sequence.fasta -p 8473 -n A -o ../data/2_derived/01KP
    python 01A_prepare_function.py --help
"""

import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def extract_motif(sequence, pos, flank=2):
    """Extract motif around a position (1‑based)."""
    start_pos = max(pos - flank - 1, 0)      # 0‑based start
    end_pos = min(pos + flank, len(sequence))
    motif = sequence[start_pos:end_pos]
    return motif, start_pos, end_pos


def hamming_distance(s1, s2):
    """Hamming distance between two equal‑length strings."""
    if len(s1) != len(s2):
        raise ValueError("Strings must have the same length")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def find_matches(sequence, motif, max_mismatch=1):
    """Find all occurrences of motif in sequence with up to max_mismatch mismatches."""
    motif_len = len(motif)
    seq_len = len(sequence)
    matches = []
    for i in range(seq_len - motif_len + 1):
        kmer = sequence[i:i+motif_len]
        dist = hamming_distance(motif, kmer)
        if dist <= max_mismatch:
            matches.append({
                'position': i + 1,          # 1‑based
                'motif': kmer,
                'distance': dist
            })
    return pd.DataFrame(matches)


def get_motif_around_position(seq, pos, left_flank, right_flank):
    """
    Extract motif with asymmetric flanks.
    Returns (motif_seq, start_1based, end_1based).
    """
    start_pos = max(pos - left_flank - 1, 0)
    end_pos = min(pos + right_flank, len(seq))
    motif = seq[start_pos:end_pos]
    return motif, start_pos + 1, end_pos


def find_approximate_repeats(seq, motif, max_mismatch):
    """Find all approximate repeats of motif in seq."""
    motif_len = len(motif)
    repeats = []
    for i in range(len(seq) - motif_len + 1):
        kmer = seq[i:i+motif_len]
        dist = hamming_distance(motif, kmer)
        if dist <= max_mismatch:
            repeats.append({
                'repeat.seq': kmer,
                'repeat.start': i + 1,
                'repeat.end': i + motif_len,
                'repeat.hamming.distance': dist
            })
    return pd.DataFrame(repeats)


def get_repeatability_table(input_seq, pos, fixed_nuc, output_folder,
                           max_flank_left=20, max_flank_right=20,
                           min_length=5, max_length=41):
    """
    Main function: compute repeatability table for a given allele and position.
    Saves result to output_folder/01KP.{pos}.{fixed_nuc}.txt and returns DataFrame.
    """
    # Create allele sequence
    seq_list = list(input_seq)
    if seq_list[pos-1] != fixed_nuc:
        seq_list[pos-1] = fixed_nuc
    seq_allele = ''.join(seq_list)

    all_results = []

    for motif_length in range(min_length, max_length + 1):
        for left_flank in range(motif_length):
            right_flank = motif_length - left_flank - 1
            if left_flank <= max_flank_left and right_flank <= max_flank_right:
                motif_seq, motif_start, motif_end = get_motif_around_position(
                    seq_allele, pos, left_flank, right_flank
                )
                max_mismatch_allowed = int(0.2 * len(motif_seq))
                repeats_df = find_approximate_repeats(seq_allele, motif_seq, max_mismatch_allowed)

                if not repeats_df.empty:
                    repeats_df['motif.seq'] = motif_seq
                    repeats_df['motif.length'] = len(motif_seq)
                    repeats_df['motif.start'] = motif_start
                    repeats_df['motif.end'] = motif_end
                    repeats_df['nuc'] = fixed_nuc
                    repeats_df['effective.length'] = repeats_df['motif.length'] - repeats_df['repeat.hamming.distance']

                    for _, row in repeats_df.iterrows():
                        all_results.append({
                            'pos': pos,
                            'nuc': row['nuc'],
                            'motif.seq': row['motif.seq'],
                            'motif.length': row['motif.length'],
                            'motif.start': row['motif.start'],
                            'motif.end': row['motif.end'],
                            'repeat.seq': row['repeat.seq'],
                            'repeat.start': row['repeat.start'],
                            'repeat.end': row['repeat.end'],
                            'repeat.hamming.distance': row['repeat.hamming.distance'],
                            'effective.length': row['effective.length']
                        })

    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df = results_df.sort_values('effective.length', ascending=False)

        out_path = Path(output_folder) / f"01KP.{pos}.{fixed_nuc}.txt"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(out_path, sep='\t', index=False)
        return results_df
    else:
        print("No repeats found")
        return pd.DataFrame()


def main():
    parser = argparse.ArgumentParser(
        description="Compute repeatability table for a single position and nucleotide",
        epilog="Example: python %(prog)s --fasta data.fasta --pos 8473 --nuc T"
    )
    parser.add_argument('-f', '--fasta', required=True, help="Path to reference FASTA file")
    parser.add_argument('-p', '--pos', type=int, required=True, help="1‑based position in the genome")
    parser.add_argument('-n', '--nuc', required=True, choices=['A','T','G','C'], help="Nucleotide to place at the position")
    parser.add_argument('-o', '--outdir', default="../data/2_derived/01KP", help="Output directory (default: ../data/2_derived/01KP)")
    parser.add_argument('--max_flank_left', type=int, default=20, help="Max left flank length (default: 20)")
    parser.add_argument('--max_flank_right', type=int, default=20, help="Max right flank length (default: 20)")
    parser.add_argument('--min_len', type=int, default=5, help="Minimum motif length (default: 5)")
    parser.add_argument('--max_len', type=int, default=41, help="Maximum motif length (default: 41)")
    parser.add_argument('--quiet', action='store_true', help="Suppress extra output (only errors)")
    args = parser.parse_args()

    # Load reference
    try:
        record = next(SeqIO.parse(args.fasta, "fasta"))
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return 1
    mtDNA_seq = str(record.seq).upper()
    if not args.quiet:
        print(f"mtDNA length: {len(mtDNA_seq)}")
        print(f"First 100 nucleotides: {mtDNA_seq[:100]}")

    # QC demonstration (similar to original notebook)
    if not args.quiet:
        flank_demo = 2
        motif_demo, start_demo, end_demo = extract_motif(mtDNA_seq, args.pos, flank_demo)
        print(f"\nMotif (positions {start_demo+1}-{end_demo}): {motif_demo}")
        print(f"Motif length: {len(motif_demo)}")
        max_mm_demo = 1
        matches_demo = find_matches(mtDNA_seq, motif_demo, max_mm_demo)
        print(f"\nFound matches: {len(matches_demo)}")
        print(f"Perfect matches: {len(matches_demo[matches_demo['distance'] == 0])}")
        print(f"Matches with 1 mismatch: {len(matches_demo[matches_demo['distance'] == 1])}")
        print("\nFirst 10 matches:")
        print(matches_demo.head(10))

    # Compute repeatability table
    results_df = get_repeatability_table(
        mtDNA_seq,
        args.pos,
        args.nuc.upper(),
        args.outdir,
        max_flank_left=args.max_flank_left,
        max_flank_right=args.max_flank_right,
        min_length=args.min_len,
        max_length=args.max_len
    )

    if not args.quiet:
        print(f"\nTotal repeats found: {len(results_df)}")
        if not results_df.empty:
            print("First 5 repeats:")
            print(results_df.head())
    return 0


if __name__ == "__main__":
    exit(main())
