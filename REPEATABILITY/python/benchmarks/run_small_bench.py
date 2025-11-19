import time
import statistics
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.analysis import get_repeatability_table


def time_func(func, repeats=2):
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        func()
        t1 = time.perf_counter()
        times.append(t1 - t0)
    return times


def run_small_synthetic(length=20000, pos=None, n_workers=2):
    print('\nSmall synthetic benchmark (length=%d, n_workers=%d)' % (length, n_workers))
    seq = ('ACGT' * (length // 4 + 1))[:length]
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
            get_repeatability_table(seq, pos, 'A', 'G', max_flank_left=20, max_flank_right=20, min_length=5, max_length=41, use_numpy=cfg['use_numpy'], use_parallel=cfg['use_parallel'], n_workers=n_workers)
        times = time_func(work, repeats=2)
        print('  times (s):', [f'{t:.3f}' for t in times], 'median:', f'{statistics.median(times):.3f}')
        results.append((name, times))
    return results


if __name__ == '__main__':
    run_small_synthetic()
