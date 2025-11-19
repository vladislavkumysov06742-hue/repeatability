import subprocess
import shutil
from pathlib import Path
import sys
import pandas as pd
import math

# Debug script to run R and Python implementations and find first mismatch
# Usage: run from repo root or anywhere; script uses paths relative to its location.

repo_root = Path(__file__).resolve().parents[2]
py_dir = repo_root / "python"

# tmp working dir inside repo to store outputs
tmp = py_dir / ".pytest_tmp" / "debug_run_tmp"
if tmp.exists():
    # clear previous files that might conflict
    for f in tmp.iterdir():
        try:
            if f.is_file():
                f.unlink()
            else:
                # skip directories
                pass
        except Exception:
            pass
else:
    tmp.mkdir(parents=True, exist_ok=True)

# Prepare small fasta
fasta = tmp / "test.fasta"
seq = ("ACGT" * 500)[:2000]
fasta.write_text(">test\n" + seq + "\n")

pos = 100
ref = "A"
alt = "G"
outpref = "my_test_repeat"

# Verify Rscript presence
if shutil.which("Rscript") is None:
    print("Rscript not found on PATH. Aborting.")
    sys.exit(2)

# Path to R functions file
r_functions = repo_root / "scripts" / "02.RefAltRepeatability.Functions.R"
if not r_functions.exists():
    print("R functions file not found:", r_functions)
    sys.exit(3)

# Create small R runner script in tmp
r_runner = tmp / "run_r_analysis.R"
r_runner.write_text(
    f"source('{r_functions.as_posix()}')\n" +
    f"analyze_mtDNA_repeats(fasta_path = '{fasta.as_posix()}', pos = {pos}, ref_nuc = '{ref}', alt_nuc = '{alt}', output_prefix = '{outpref}')\n"
)

print('Running Rscript...')
subprocess.run(["Rscript", r_runner.as_posix()], check=True, cwd=str(tmp))
print('R finished.')

# Run Python implementation
sys.path.insert(0, str(py_dir))
from src.analysis import analyze_mtDNA_repeats
print('Running Python analyze_mtDNA_repeats...')
analyze_mtDNA_repeats(str(fasta), pos, ref, alt, output_path=str(tmp), output_prefix=outpref)
print('Python finished.')

# Paths to outputs
rA = tmp / f"01KP.{pos}.{ref}.txt"
rG = tmp / f"01KP.{pos}.{alt}.txt"
py_dir_out = tmp / f"pos{pos}_{ref}to{alt}"
pyA = py_dir_out / f"01KP.{pos}.{ref}.txt"
pyG = py_dir_out / f"01KP.{pos}.{alt}.txt"

print('\nFiles produced:')
for p in [rA, rG, pyA, pyG]:
    print(p, 'exists' if p.exists() else 'MISSING')


def load_tab(p):
    if not p.exists():
        print('Missing', p)
        return pd.DataFrame()
    try:
        return pd.read_csv(p, sep="\t", dtype=str)
    except Exception as e:
        print('Failed to read', p, e)
        return pd.DataFrame()

# Compare for each allele
for allele, r_file, py_file in ((ref, rA, pyA), (alt, rG, pyG)):
    print(f"\nComparing allele {allele}:")
    r_df = load_tab(r_file)
    py_df = load_tab(py_file)
    print('R rows:', len(r_df), 'PY rows:', len(py_df))
    print('R cols:', list(r_df.columns))
    print('PY cols:', list(py_df.columns))
    if r_df.empty and py_df.empty:
        print('Both empty, skipping')
        continue

    # Align columns
    all_cols = sorted(set(r_df.columns).union(set(py_df.columns)))
    for c in all_cols:
        if c not in r_df.columns:
            r_df[c] = ""
        if c not in py_df.columns:
            py_df[c] = ""

    r_s = r_df[all_cols].astype(str).sort_values(by=all_cols).reset_index(drop=True)
    py_s = py_df[all_cols].astype(str).sort_values(by=all_cols).reset_index(drop=True)

    # Quick equality check
    try:
        pd.testing.assert_frame_equal(r_s, py_s, check_like=True, check_dtype=False)
        print('DataFrames match exactly for allele', allele)
        continue
    except AssertionError as e:
        print('DataFrames differ for allele', allele)

    # Find first mismatch
    # pad shorter with empty rows
    max_rows = max(len(r_s), len(py_s))
    if len(r_s) < max_rows:
        r_s = pd.concat([r_s, pd.DataFrame([[''] * len(all_cols)] * (max_rows - len(r_s),), columns=all_cols)], ignore_index=True)
    if len(py_s) < max_rows:
        py_s = pd.concat([py_s, pd.DataFrame([[''] * len(all_cols)] * (max_rows - len(py_s),), columns=all_cols)], ignore_index=True)

    mismatch = (r_s.values != py_s.values)
    any_mismatch = mismatch.any()
    if not any_mismatch:
        print('No cell-wise mismatches after padding — possible dtype/index differences')
        continue
    # find first True
    idx = int(mismatch.argmax())
    row = idx // mismatch.shape[1]
    col = idx % mismatch.shape[1]
    print('First mismatch at row', row, 'column', all_cols[col])
    print('R value:', repr(r_s.iloc[row, col]))
    print('PY value:', repr(py_s.iloc[row, col]))

    # Show context
    start = max(0, row-3)
    end = min(max_rows, row+4)
    print('\nR rows around mismatch:')
    print(r_s.iloc[start:end].to_string(index=True))
    print('\nPY rows around mismatch:')
    print(py_s.iloc[start:end].to_string(index=True))

print('\nDebug run complete.')
