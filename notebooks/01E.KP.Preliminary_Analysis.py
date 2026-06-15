#!/usr/bin/env python3
"""
01E.KP.Preliminary_Analysis.py

Robust repeatability metric extraction for specified positions.
For each position and nucleotide, reads the corresponding 01KP.{pos}.{nuc}.txt file,
applies filtering (major arc, exclude self-match, perfect repeats only, pos inside motif),
and returns the longest perfect repeat length (or 0 if none).

Author: Vladislav Gadzhiev
Date: 13/12/2025

Usage examples:
    python 01E_preliminary_analysis.py --input_dir ../data/2_derived/01KP --positions 8251 8473 12705
    python 01E_preliminary_analysis.py -i ../data/2_derived/01KP -p 8251 8472 8473 -o results.csv --nucleotides A T G C
    python 01E_preliminary_analysis.py --help
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def safe_read_repeatability_file(pos, nuc, results_dir):
    """
    Safely read a repeatability file. If the file does not exist or is invalid,
    returns a minimal DataFrame with NaN values for that combination.
    """
    file_path = Path(results_dir) / f"01KP.{pos}.{nuc}.txt"
    try:
        df = pd.read_csv(file_path, sep='\t')
        # Add position and nucleotide columns for convenience
        df['pos'] = pos
        df['nuc'] = nuc
        return df
    except FileNotFoundError:
        print(f"  File not found, skipping: {file_path.name}")
        # Return a minimal row with NaNs so that we can still have an entry
        empty_data = {
            'pos': pos,
            'nuc': nuc,
            'motif.length': np.nan,
            'effective.length': np.nan,
            'repeat.hamming.distance': np.nan,
            'motif.start': np.nan,
            'motif.end': np.nan,
            'repeat.start': np.nan,
            'repeat.end': np.nan
        }
        return pd.DataFrame([empty_data])
    except Exception as e:
        print(f"  Error reading {file_path.name}: {e}")
        return pd.DataFrame()  # empty DataFrame, will be skipped


def calculate_robust_repeatability_metric(results_dir, positions_of_interest,
                                          nucleotides=None,
                                          major_arc_start=5798,
                                          major_arc_end=16568,
                                          output_csv=None):
    """
    Compute the longest perfect repeat length for given positions and nucleotides.
    Filters:
        - motif != repeat (exact start/end match)
        - both motif and repeat lie within the major arc
        - repeat.hamming.distance == 0 (perfect match)
        - pos strictly inside motif (motif.start < pos < motif.end)
    Returns a DataFrame with columns: position, nucleotide, longest_perfect_repeat_length.
    If output_csv is provided, saves the result.
    """
    if nucleotides is None:
        nucleotides = ['A', 'T', 'G', 'C']

    results = []
    total = len(positions_of_interest) * len(nucleotides)
    processed = 0

    for pos in positions_of_interest:
        for nuc in nucleotides:
            df = safe_read_repeatability_file(pos, nuc, results_dir)
            processed += 1
            if processed % 100 == 0:
                print(f"  Progress: {processed}/{total}")

            if df.empty:
                continue

            # Apply filters exactly as in the original notebook
            filtered = df[
                ~((df['motif.start'] == df['repeat.start']) & (df['motif.end'] == df['repeat.end'])) &
                (df['motif.start'] >= major_arc_start) & (df['motif.end'] <= major_arc_end) &
                (df['repeat.start'] >= major_arc_start) & (df['repeat.end'] <= major_arc_end) &
                (df['repeat.hamming.distance'] == 0) &
                (df['motif.start'] < df['pos']) & (df['pos'] < df['motif.end'])
            ]

            if not filtered.empty:
                longest_length = filtered['motif.length'].max()
            else:
                longest_length = 0  # no perfect repeat meeting criteria

            results.append({
                'position': pos,
                'nucleotide': nuc,
                'longest_perfect_repeat_length': longest_length
            })

    result_df = pd.DataFrame(results)

    if output_csv:
        output_path = Path(output_csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output_path, index=False)
        print(f"Results saved to {output_path}")

    return result_df


def main():
    parser = argparse.ArgumentParser(
        description="Compute longest perfect repeat length for specified positions.",
        epilog="Example: %(prog)s --input_dir data/2_derived/01KP --positions 8251 8473 12705"
    )
    parser.add_argument('-i', '--input_dir', required=True,
                        help="Directory containing 01KP.*.txt files")
    parser.add_argument('-p', '--positions', type=int, nargs='+', required=True,
                        help="List of 1‑based positions to analyze (e.g., 8251 8473 12705)")
    parser.add_argument('-n', '--nucleotides', nargs='+', default=['A','T','G','C'],
                        choices=['A','T','G','C'],
                        help="Nucleotides to test (default: A T G C)")
    parser.add_argument('--major_start', type=int, default=5798,
                        help="Start of major arc (default: 5798)")
    parser.add_argument('--major_end', type=int, default=16568,
                        help="End of major arc (default: 16568)")
    parser.add_argument('-o', '--output', default='repeatability_metrics.csv',
                        help="Output CSV file (default: repeatability_metrics.csv)")
    args = parser.parse_args()

    # Convert input_dir to Path object
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"Error: input directory {input_dir} does not exist.")
        return 1

    print(f"Analyzing {len(args.positions)} positions with nucleotides {args.nucleotides}")
    df = calculate_robust_repeatability_metric(
        results_dir=input_dir,
        positions_of_interest=args.positions,
        nucleotides=args.nucleotides,
        major_arc_start=args.major_start,
        major_arc_end=args.major_end,
        output_csv=args.output
    )

    # Print summary
    if not df.empty:
        print("\nSummary:")
        print(f"  Total entries: {len(df)}")
        non_zero = df[df['longest_perfect_repeat_length'] > 0]
        print(f"  Positions with at least one perfect repeat: {non_zero['position'].nunique()}")
        print(f"  Mean longest length (non-zero): {non_zero['longest_perfect_repeat_length'].mean():.2f}")
        print("\nFirst few rows:")
        print(df.head(10))
    return 0


if __name__ == "__main__":
    exit(main())
