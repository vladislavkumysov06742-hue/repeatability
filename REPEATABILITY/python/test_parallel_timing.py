#!/usr/bin/env python3
import time
import subprocess
import sys

# Test parallel (2 workers) on 3 positions
start = time.time()
result = subprocess.run(
    [sys.executable, 'run_all_variants.py', '--start', '8006', '--end', '8008', '--use-numpy', '--max-workers', '2'],
    capture_output=True,
    text=True
)
parallel_time = time.time() - start

# Test sequential (1 worker) on 3 positions
start = time.time()
result = subprocess.run(
    [sys.executable, 'run_all_variants.py', '--start', '8009', '--end', '8011', '--use-numpy', '--max-workers', '1'],
    capture_output=True,
    text=True
)
sequential_time = time.time() - start

print("\n" + "="*60)
print(f"Parallel (2 workers, 3 positions):   {parallel_time:.2f} sec")
print(f"Sequential (1 worker, 3 positions):  {sequential_time:.2f} sec")
print(f"Speedup: {sequential_time/parallel_time:.2f}x")
print("="*60)
