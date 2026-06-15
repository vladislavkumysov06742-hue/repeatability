#!/usr/bin/env python3
"""
01B.KP.RunTheFunctionOnWholeMtDna.py

Calculate repeatability for Ref and Alt alleles for all positions (or a range)
within human mtDNA, using parallel processing.

Author: Vladislav Gadzhiev
Date: 11/12/2025

Usage examples:
    python 01B_run_whole_mtDNA.py --fasta ../data/1_raw/sequence.fasta
    python 01B_run_whole_mtDNA.py -f data.fasta -o results --start 1 --end 1000 --workers 8
    python 01B_run_whole_mtDNA.py --help
"""

import argparse
import time
import numpy as np
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import warnings
warnings.filterwarnings('ignore')


# ---------- Core repeatability functions (copied from 01A for self‑containment) ----------
def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strings must have the same length")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def get_motif_around_position(seq, pos, left_flank, right_flank):
    start_pos = max(pos - left_flank - 1, 0)
    end_pos = min(pos + right_flank, len(seq))
    motif = seq[start_pos:end_pos]
    return motif, start_pos + 1, end_pos


def find_approximate_repeats(seq, motif, max_mismatch):
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
    """Compute and save repeatability table for a single allele."""
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
        return pd.DataFrame()


# ---------- Parallel processing wrapper ----------
def process_single(args):
    """Wrapper for get_repeatability_table to be used with ProcessPoolExecutor."""
    seq, pos, nuc, outdir = args
    df = get_repeatability_table(seq, pos, nuc, outdir)
    return (pos, nuc, len(df))


def main():
    parser = argparse.ArgumentParser(
        description="Run repeatability calculation for a range of positions (default: whole genome).",
        epilog="Example: python %(prog)s --fasta mtDNA.fasta --outdir results --workers 16"
    )
    parser.add_argument('-f', '--fasta', required=True, help="Path to reference FASTA file")
    parser.add_argument('-o', '--outdir', default="../data/2_derived/01KP", help="Output directory (default: ../data/2_derived/01KP)")
    parser.add_argument('--start', type=int, default=1, help="First position to process (1‑based, default: 1)")
    parser.add_argument('--end', type=int, help="Last position to process (default: genome length)")
    parser.add_argument('--workers', type=int, default=mp.cpu_count(), help=f"Number of parallel workers (default: CPU count = {mp.cpu_count()})")
    parser.add_argument('--batch_size', type=int, default=100, help="Number of positions per batch (for progress reporting, default: 100)")
    parser.add_argument('--nucleotides', nargs='+', default=['A','T','G','C'], choices=['A','T','G','C'], help="Nucleotides to test (default: A T G C)")
    parser.add_argument('--quiet', action='store_true', help="Suppress detailed progress output")
    args = parser.parse_args()

    # Load reference sequence
    record = next(SeqIO.parse(args.fasta, "fasta"))
    seq = str(record.seq).upper()
    genome_len = len(seq)
    end_pos = args.end if args.end is not None else genome_len

    if args.start < 1 or end_pos > genome_len or args.start > end_pos:
        print(f"Error: invalid position range ({args.start}-{end_pos}). Genome length = {genome_len}")
        return 1

    total_positions = end_pos - args.start + 1
    total_tasks = total_positions * len(args.nucleotides)
    if not args.quiet:
        print(f"Genome length: {genome_len}")
        print(f"Processing positions {args.start} to {end_pos} ({total_positions} positions)")
        print(f"Nucleotides: {args.nucleotides}")
        print(f"Total tasks: {total_tasks}")
        print(f"Workers: {args.workers}")
        print(f"Output directory: {args.outdir}")

    # Create output directory
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # Prepare all tasks
    tasks = [(seq, pos, nuc, args.outdir) for pos in range(args.start, end_pos+1) for nuc in args.nucleotides]

    # Process in batches for progress reporting (optional)
    if args.batch_size <= 0:
        args.batch_size = total_positions

    start_time = time.time()
    processed = 0

    # Split tasks into batches of size batch_size * len(nucleotides)
    batch_size_tasks = args.batch_size * len(args.nucleotides)
    num_batches = (len(tasks) + batch_size_tasks - 1) // batch_size_tasks

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for batch_idx in range(num_batches):
            batch_tasks = tasks[batch_idx * batch_size_tasks : (batch_idx+1) * batch_size_tasks]
            if not batch_tasks:
                break
            # Submit batch
            futures = [executor.submit(process_single, t) for t in batch_tasks]
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    print(f"Error in task: {e}")
                processed += 1
                if not args.quiet and processed % 100 == 0:
                    elapsed = time.time() - start_time
                    rate = processed / elapsed if elapsed > 0 else 0
                    print(f"Progress: {processed}/{total_tasks} ({100*processed/total_tasks:.1f}%) "
                          f"- {rate:.1f} tasks/sec, elapsed {elapsed:.1f}s")

    elapsed = time.time() - start_time
    if not args.quiet:
        print(f"\n{'='*50}")
        print(f"ALL TASKS COMPLETED")
        print(f"Total time: {elapsed:.1f} sec")
        print(f"Tasks processed: {processed}")
        print(f"Average time per task: {elapsed/processed:.3f} sec" if processed else "")
    return 0


if __name__ == "__main__":
    exit(main())
