import time
import statistics
from pathlib import Path
import argparse
import sys

# ensure package import
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.analysis import analyze_mtDNA_repeats, get_repeatability_table


def time_func(func, repeats=3):
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        func()
        t1 = time.perf_counter()
        times.append(t1 - t0)
    return times


def run_bench_real(fasta_path, pos=8473, ref='T', alt='C'):
    print('Benchmark on real mtDNA fasta:', fasta_path)
    configs = [
        {'name': 'pure-python', 'use_numpy': False, 'use_parallel': False},
        {'name': 'numpy', 'use_numpy': True, 'use_parallel': False},
        {'name': 'parallel-pure', 'use_numpy': False, 'use_parallel': True},
        {'name': 'parallel-numpy', 'use_numpy': True, 'use_parallel': True},
    ]
    results = []
    for cfg in configs:
        name = cfg['name']
        print('\nRunning config:', name)
        def work():
            analyze_mtDNA_repeats(str(fasta_path), pos, ref, alt, output_path='.', use_numpy=cfg['use_numpy'], use_parallel=cfg['use_parallel'])
        times = time_func(work, repeats=3)
        print('  times (s):', [f'{t:.3f}' for t in times], 'median:', f'{statistics.median(times):.3f}')
        results.append((name, times))
    return results


def run_bench_synthetic(length=100000, motif='ATGCA', pos=None):
    print('\nBenchmark on synthetic sequence (length=%d)' % length)
    # build a synthetic sequence with some repeats
    seq = ('ACGT' * (length // 4 + 1))[:length]
    # create a fasta-like file for analyze_mtDNA_repeats which expects a real fasta path; instead we will use get_repeatability_table directly
    # We will time get_repeatability_table with different settings but using the same in-memory sequence
    pos = pos or (length // 2)
    configs = [
        {'name': 'pure-python', 'use_numpy': False, 'use_parallel': False},
        {'name': 'numpy', 'use_numpy': True, 'use_parallel': False},
        {'name': 'parallel-pure', 'use_numpy': False, 'use_parallel': True},
        {'name': 'parallel-numpy', 'use_numpy': True, 'use_parallel': True},
    ]
    results = []
    for cfg in configs:
        name = cfg['name']
        print('\nRunning config:', name)
        def work():
            # call get_repeatability_table directly to avoid fasta I/O
            get_repeatability_table(seq, pos, 'A', 'G', max_flank_left=20, max_flank_right=20, min_length=5, max_length=41, use_numpy=cfg['use_numpy'], use_parallel=cfg['use_parallel'])
        times = time_func(work, repeats=3)
        print('  times (s):', [f'{t:.3f}' for t in times], 'median:', f'{statistics.median(times):.3f}')
        results.append((name, times))
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Benchmark repeat search configurations')
    parser.add_argument('--fasta', type=str, default='../data/1_raw/Homo_sapients.mtDNA.fasta')
    parser.add_argument('--pos', type=int, default=8473)
    parser.add_argument('--synthetic-length', type=int, default=50000)
    args = parser.parse_args()

    fasta = Path(args.fasta).resolve()
    if fasta.exists():
        run_bench_real(fasta, pos=args.pos)
    else:
        print('Real fasta not found at', fasta)

    run_bench_synthetic(length=args.synthetic_length)
