#!/usr/bin/env python3
import time
import subprocess
import sys

positions = [8251, 8472, 8473]
times = {}

print("Testing real biological positions with ALL optimizations...")
print("="*60)

for pos in positions:
    start = time.time()
    result = subprocess.run(
        [sys.executable, 'run_all_variants.py', '--start', str(pos), '--end', str(pos), '--use-numpy'],
        capture_output=True,
        text=True
    )
    elapsed = time.time() - start
    times[pos] = elapsed
    
    status = "✅" if result.returncode == 0 else "❌"
    print(f"{status} pos{pos}: {elapsed:.2f} sec")

total = sum(times.values())
avg_per_pos = total / len(positions)

print("="*60)
print(f"Total (3 real positions): {total:.2f} sec")
print(f"Average per position: {avg_per_pos:.2f} sec")
print(f"\nComparison to baseline:")
print(f"  Without optimizations: ~65 sec per position")
print(f"  With all optimizations: {avg_per_pos:.2f} sec per position")
print(f"  Speedup: {65/avg_per_pos:.1f}x")
