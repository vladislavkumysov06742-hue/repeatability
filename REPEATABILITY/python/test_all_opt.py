#!/usr/bin/env python3
import time
import subprocess
import sys

# Test with all optimizations (3 positions)
print("Testing with ALL optimizations (numpy + batch DF + loaded seq)...")
start = time.time()
result = subprocess.run(
    [sys.executable, 'run_all_variants.py', '--start', '8015', '--end', '8017', '--use-numpy'],
    capture_output=True,
    text=True
)
elapsed_all_opt = time.time() - start

# Print result
if result.returncode == 0:
    print(f"\n✅ All optimizations (3 positions, 9 analyses): {elapsed_all_opt:.2f} sec")
    print(f"   Average per position: {elapsed_all_opt/3:.2f} sec")
    print(f"   Average per alt: {elapsed_all_opt/9:.2f} sec")
else:
    print(f"❌ Error: {result.stderr}")
    
# Compare with baseline (no optimizations would be ~65 sec for 3 positions; 
# with numpy alone ~21 sec each = 63 sec; with all opts should be ~38-40 sec)
print("\n" + "="*60)
print("Expected baselines:")
print("  No optimizations: ~65 sec (for 3 pos)")
print("  NumPy only: ~63 sec")
print("  All opts (numpy + batch + seq): ~38-42 sec (expected)")
print("="*60)
