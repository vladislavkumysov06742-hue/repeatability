#!/usr/bin/env python3
"""Run analyze_mtDNA_repeats for every position and all alternative alleles.

This is the runner we used interactively. It supports --use-numpy, --precompute-ref-only
and --max-workers to control parallelism.
"""
from pathlib import Path
import argparse
import sys


def read_fasta_sequence(fasta_path: Path) -> str:
    seq = []
    with open(fasta_path, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                continue
            seq.append(line.strip())
    return ''.join(seq)


def main():
    parser = argparse.ArgumentParser(description='Run repeat analysis for all alt alleles at all positions')
    parser.add_argument('--fasta', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / '1_raw' / 'Homo_sapients.mtDNA.fasta',
                        help='Reference FASTA (default: data/1_raw/Homo_sapients.mtDNA.fasta)')
    parser.add_argument('--outdir', type=Path, default=Path(__file__).resolve().parents[1] / 'data' / '2_derived',
                        help='Output base directory for per-pos results')
    parser.add_argument('--start', type=int, default=1, help='Start position (1-based)')
    parser.add_argument('--end', type=int, default=None, help='End position (1-based, inclusive). If omitted uses full length')
    parser.add_argument('--skip-existing', action='store_true', help='Skip analyses if output folder exists')
    parser.add_argument('--dry-run', action='store_true', help='Do not run analysis; just print planned actions')
    parser.add_argument('--max', type=int, default=None, help='Limit to first N positions (for quick tests)')
    parser.add_argument('--use-numpy', action='store_true', default=False, help='Use NumPy-vectorized search (faster)')
    parser.add_argument('--precompute-ref-only', action='store_true', help='Only precompute reference for positions')
    parser.add_argument('--max-workers', type=int, default=1, help='Max workers for parallel position processing (default=1=sequential)')
    args = parser.parse_args()

    fasta = args.fasta
    if not fasta.exists():
        print('Fasta not found:', fasta, file=sys.stderr)
        sys.exit(2)

    seq = read_fasta_sequence(fasta)
    length = len(seq)
    end = args.end if args.end is not None else length
    if args.start < 1 or end > length or args.start > end:
        print(f'Invalid start/end: start={args.start}, end={end}, sequence length={length}', file=sys.stderr)
        sys.exit(2)

    # Import the analysis function lazily so the script can be inspected without heavy deps
    try:
        from src.analysis import analyze_mtDNA_repeats
    except Exception as e:
        print('Failed to import analyze_mtDNA_repeats from src.analysis:', e, file=sys.stderr)
        if not args.dry_run:
            sys.exit(3)

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    cache_dir_arg = str(outdir / 'cache')

    positions = list(range(args.start, end + 1))
    if args.max is not None:
        import itertools
        positions = list(itertools.islice(positions, args.max))

    base_order = ['T', 'G', 'C', 'A']

    if args.precompute_ref_only:
        print('Precomputing references only...')
        for pos in positions:
            ref = seq[pos - 1].upper()
            if ref not in {'A', 'C', 'G', 'T'}:
                continue
            print(f'PRECOMPUTE REF pos {pos} {ref}')
            if args.dry_run:
                continue
            try:
                analyze_mtDNA_repeats(str(fasta), pos, ref, ref, output_path=str(outdir), output_prefix=f'mt_repeat_pos{pos}', compute_ref=True, cache_dir=cache_dir_arg, use_numpy=args.use_numpy, seq=seq)
            except Exception as e:
                print(f'ERROR precomputing ref for pos {pos}:', e, file=sys.stderr)
        print('Precompute done. Run without --precompute-ref-only to process alternatives.')
        return

    total = [0]

    def process_position(pos):
        ref = seq[pos - 1].upper()
        if ref not in {'A', 'C', 'G', 'T'}:
            return
        alt_list = [n for n in base_order if n != ref]
        for i, alt in enumerate(alt_list):
            total[0] += 1
            folder = outdir / f'pos{pos}_{ref}to{alt}'
            if args.skip_existing and folder.exists():
                print(f'SKIP pos {pos} {ref}->{alt} (exists)')
                continue
            print(f'RUN pos {pos} {ref}->{alt} -> {folder}')
            if args.dry_run:
                continue
            compute_ref = (i == 0)
            try:
                analyze_mtDNA_repeats(str(fasta), pos, ref, alt, output_path=str(outdir), output_prefix=f'mt_repeat_pos{pos}', compute_ref=compute_ref, cache_dir=cache_dir_arg, use_numpy=args.use_numpy, seq=seq)
            except Exception as e:
                print(f'ERROR running pos {pos} {ref}->{alt}:', e, file=sys.stderr)

    if args.max_workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        print(f'Processing positions with {args.max_workers} workers...')
        with ProcessPoolExecutor(max_workers=args.max_workers) as exe:
            list(exe.map(process_position, positions))
    else:
        for pos in positions:
            process_position(pos)

    print(f'Done. Planned/attempted analyses: {total[0]}')


if __name__ == '__main__':
    main()
