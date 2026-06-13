#!/usr/bin/env python3
"""Run analyze_mtDNA_repeats for every position and all alternative alleles.

Usage examples:
  python run_all_variants.py --start 1 --end 100 --skip-existing
  python run_all_variants.py --dry-run --start 8000 --end 9000

The script calls `retrieve_all_repeats_covering_given_position_and_nucleotide_in_sequence`
for each position and nucleotide (A, T, G, C) and saves results to output files.
"""

import argparse
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from typing import List, Tuple, Optional
import time
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent / "functions"))

try:
    from RetrieveAllRepeatsCoveringGivenPositionAndNucleotideInSequence import (
        retrieve_all_repeats_covering_given_position_and_nucleotide_in_sequence,
    )
except ImportError:
    print("Error: Could not import repeat analysis functions.", file=sys.stderr)
    print(
        "Make sure the function file is in the 'functions' directory.", file=sys.stderr
    )
    sys.exit(1)


def read_fasta_sequence(fasta_path: Path) -> str:
    """Read DNA sequence from FASTA file."""
    seq = []
    with open(fasta_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


def process_single_position(
    seq: str,
    pos: int,
    nuc: str,
    output_folder: Path,
    skip_existing: bool = False,
    dry_run: bool = False,
) -> Tuple[int, str, bool, str]:
    """
    Process a single position-nucleotide combination.

    Returns: (position, nucleotide, success, message)
    """
    output_file = output_folder / f"01KP.{pos}.{nuc}.txt"

    if skip_existing and output_file.exists():
        return (pos, nuc, True, f"SKIP (exists): {output_file.name}")

    if dry_run:
        return (pos, nuc, True, f"DRY-RUN: Would process {pos} {nuc}")

    try:
        results = (
            retrieve_all_repeats_covering_given_position_and_nucleotide_in_sequence(
                seq=seq,
                pos=pos,
                fixed_nuc=nuc,
                output_folder=str(output_folder),
                max_flank_left=20,
                max_flank_right=20,
                min_length=5,
                max_length=41,
                mismatch_fraction=0.2,
            )
        )

        if len(results) > 0:
            return (pos, nuc, True, f"SUCCESS: {len(results)} repeats found")
        else:
            return (pos, nuc, True, f"SUCCESS: No repeats found (file created)")

    except Exception as e:
        return (pos, nuc, False, f"ERROR: {str(e)}")


def create_tasks(
    seq: str,
    positions: List[int],
    nucleotides: List[str],
    output_folder: Path,
    skip_existing: bool,
    dry_run: bool,
) -> List[dict]:
    """Create task list for processing."""
    tasks = []
    for pos in positions:
        for nuc in nucleotides:
            tasks.append(
                {
                    "seq": seq,
                    "pos": pos,
                    "nuc": nuc,
                    "output_folder": output_folder,
                    "skip_existing": skip_existing,
                    "dry_run": dry_run,
                }
            )
    return tasks


def worker_wrapper(args: dict) -> Tuple[int, str, bool, str]:
    """Wrapper function for parallel processing."""
    return process_single_position(**args)


def main():
    parser = argparse.ArgumentParser(
        description="Run repeat analysis for all alt alleles at all positions in mtDNA"
    )

    parser.add_argument(
        "--fasta",
        type=Path,
        default=Path(__file__).resolve().parent
        / "data"
        / "1_raw"
        / "Homo_sapients.mtDNA.fasta",
        help="Reference FASTA (default: data/1_raw/Homo_sapients.mtDNA.fasta)",
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parent / "data" / "2_derived" / "01KP",
        help="Output base directory for results",
    )

    parser.add_argument(
        "--start", type=int, default=1, help="Start position (1-based, inclusive)"
    )

    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help="End position (1-based, inclusive). If omitted uses full sequence length",
    )

    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip analyses if output file exists",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Do not run analysis; just print planned actions",
    )

    parser.add_argument(
        "--max-positions",
        type=int,
        default=None,
        dest="max",
        help="Limit to first N positions (for quick tests)",
    )

    parser.add_argument(
        "--max-workers",
        type=int,
        default=1,
        help="Max workers for parallel processing (default=1=sequential)",
    )

    parser.add_argument(
        "--nucleotides",
        type=str,
        default="ATGC",
        help="Nucleotides to process (default: ATGC for all)",
    )

    parser.add_argument(
        "--progress", action="store_true", help="Show progress bar (requires tqdm)"
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=1000,
        help="Number of positions to process in memory before saving (for optimization)",
    )

    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Log file to save processing results",
    )

    args = parser.parse_args()

    if not args.fasta.exists():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(2)

    print(f"Reading FASTA from: {args.fasta}")
    seq = read_fasta_sequence(args.fasta)
    seq_length = len(seq)

    end_pos = args.end if args.end is not None else seq_length

    if args.start < 1 or args.start > seq_length:
        print(
            f"Error: Start position must be between 1 and {seq_length}", file=sys.stderr
        )
        sys.exit(2)

    if end_pos < args.start or end_pos > seq_length:
        print(
            f"Error: End position must be between {args.start} and {seq_length}",
            file=sys.stderr,
        )
        sys.exit(2)

    args.outdir.mkdir(parents=True, exist_ok=True)

    nucleotides = list(args.nucleotides.upper())
    valid_nucleotides = {"A", "T", "G", "C"}

    for nuc in nucleotides:
        if nuc not in valid_nucleotides:
            print(
                f'Error: Invalid nucleotide "{nuc}". Must be one of: A, T, G, C',
                file=sys.stderr,
            )
            sys.exit(2)

    positions = list(range(args.start, end_pos + 1))
    if args.max is not None:
        positions = positions[: args.max]

    total_tasks = len(positions) * len(nucleotides)

    print(f"\n{'=' * 60}")
    print(f"Sequence length: {seq_length}")
    print(f"Processing range: {args.start}-{end_pos} ({len(positions)} positions)")
    print(f"Nucleotides: {', '.join(nucleotides)}")
    print(f"Total tasks: {total_tasks}")
    print(f"Output directory: {args.outdir}")
    print(f"Skip existing: {args.skip_existing}")
    print(f"Dry run: {args.dry_run}")
    print(f"Max workers: {args.max_workers}")
    print(f"{'=' * 60}\n")

    if args.dry_run:
        print("DRY RUN MODE - No actual processing will be done")

    tasks = create_tasks(
        seq=seq,
        positions=positions,
        nucleotides=nucleotides,
        output_folder=args.outdir,
        skip_existing=args.skip_existing,
        dry_run=args.dry_run,
    )

    start_time = time.time()
    results = []

    if args.max_workers > 1 and len(tasks) > 1:
        print(f"Starting parallel processing with {args.max_workers} workers...")

        with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
            future_to_task = {
                executor.submit(worker_wrapper, task): task for task in tasks
            }

            completed = 0
            for future in as_completed(future_to_task):
                completed += 1
                result = future.result()
                results.append(result)

                percent = (completed / total_tasks) * 100
                elapsed = time.time() - start_time

                if completed % 10 == 0 or completed == total_tasks:
                    print(
                        f"\rProgress: {completed}/{total_tasks} ({percent:.1f}%) | "
                        f"Elapsed: {elapsed:.1f}s | "
                        f"ETA: {(elapsed / completed) * (total_tasks - completed):.1f}s"
                        if completed > 0
                        else "Progress: ...",
                        end="",
                        flush=True,
                    )

    else:
        print("Starting sequential processing...")

        for i, task in enumerate(tasks, 1):
            result = worker_wrapper(task)
            results.append(result)

            percent = (i / total_tasks) * 100
            elapsed = time.time() - start_time

            if i % 10 == 0 or i == total_tasks:
                print(
                    f"\rProgress: {i}/{total_tasks} ({percent:.1f}%) | "
                    f"Elapsed: {elapsed:.1f}s | "
                    f"ETA: {(elapsed / i) * (total_tasks - i):.1f}s"
                    if i > 0
                    else "Progress: ...",
                    end="",
                    flush=True,
                )

    total_time = time.time() - start_time

    print(f"\n\n{'=' * 60}")
    print(f"PROCESSING COMPLETE")
    print(f"{'=' * 60}")
    print(f"Total time: {total_time:.2f} seconds")
    print(f"Average time per task: {total_time / len(tasks):.3f} seconds")

    success_count = sum(1 for _, _, success, _ in results if success)
    error_count = len(results) - success_count

    print(f"\nResults summary:")
    print(f"  Total tasks: {len(tasks)}")
    print(f"  Successful: {success_count}")
    print(f"  Errors: {error_count}")

    if error_count > 0:
        print(f"\nError details:")
        for pos, nuc, success, message in results:
            if not success:
                print(f"  Position {pos}, nucleotide {nuc}: {message}")

    if args.log_file:
        log_df = pd.DataFrame(
            results, columns=["position", "nucleotide", "success", "message"]
        )
        log_df.to_csv(args.log_file, sep="\t", index=False)
        print(f"\nLog saved to: {args.log_file}")

    summary_file = args.outdir / "processing_summary.txt"
    with open(summary_file, "w") as f:
        f.write(f"Processing Summary\n")
        f.write(f"{'=' * 40}\n")
        f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"FASTA: {args.fasta}\n")
        f.write(f"Sequence length: {seq_length}\n")
        f.write(f"Positions processed: {args.start}-{end_pos}\n")
        f.write(f"Nucleotides: {', '.join(nucleotides)}\n")
        f.write(f"Total tasks: {len(tasks)}\n")
        f.write(f"Successful: {success_count}\n")
        f.write(f"Errors: {error_count}\n")
        f.write(f"Total time: {total_time:.2f} seconds\n")
        f.write(f"Average time per task: {total_time / len(tasks):.3f} seconds\n")
        f.write(f"\nParameters:\n")
        f.write(f"  Skip existing: {args.skip_existing}\n")
        f.write(f"  Dry run: {args.dry_run}\n")
        f.write(f"  Max workers: {args.max_workers}\n")

    print(f"\nSummary saved to: {summary_file}")


if __name__ == "__main__":
    main()
