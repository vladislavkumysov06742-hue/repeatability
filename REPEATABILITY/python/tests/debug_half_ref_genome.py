import subprocess
import shutil
from pathlib import Path
import sys
import pandas as pd

# Script: run R and Python analysis on half of the reference mtDNA sequence
repo_root = Path(__file__).resolve().parents[2]
py_dir = repo_root / "python"

# tmp working dir
tmp = py_dir / ".pytest_tmp" / "debug_half_tmp"
if tmp.exists():
    for f in tmp.iterdir():
        try:
            if f.is_file():
                f.unlink()
            else:
                pass
        except Exception:
            pass
else:
    tmp.mkdir(parents=True, exist_ok=True)

# Build half-reference FASTA
mt_fa = repo_root / "data" / "1_raw" / "Homo_sapients.mtDNA.fasta"
if not mt_fa.exists():
    print("Reference fasta not found:", mt_fa)
    sys.exit(3)

seq_lines = []
with mt_fa.open('r', encoding='utf-8') as fh:
    for line in fh:
        if line.startswith('>'):
            continue
        seq_lines.append(line.strip())
full_seq = ''.join(seq_lines)
half_len = len(full_seq) // 2
half_seq = full_seq[:half_len]
print(f"Using half reference (len={half_len}) from {mt_fa}")

fasta = tmp / "test_half.fasta"
fasta.write_text(">test_half\n" + half_seq + "\n")

pos = 100
ref = full_seq[pos-1] if pos-1 < len(full_seq) else 'A'
alt = 'G' if ref != 'G' else 'A'
outpref = "half_ref_test"

# check Rscript
if shutil.which("Rscript") is None:
    print("Rscript not found on PATH. Aborting.")
    sys.exit(2)

# prepare R runner
r_functions = repo_root / "scripts" / "02.RefAltRepeatability.Functions.R"
if not r_functions.exists():
    print("R functions file not found:", r_functions)
    sys.exit(3)

r_runner = tmp / "run_r_analysis_half.R"
r_runner.write_text(
    f"source('{r_functions.as_posix()}')\n" +
    f"analyze_mtDNA_repeats(fasta_path = '{fasta.as_posix()}', pos = {pos}, ref_nuc = '{ref}', alt_nuc = '{alt}', output_prefix = '{outpref}')\n"
)

print('Running Rscript on half-genome...')
subprocess.run(["Rscript", r_runner.as_posix()], check=True, cwd=str(tmp))
print('R finished.')

# Run Python analyze (with debug)
sys.path.insert(0, str(py_dir))
from src.analysis import analyze_mtDNA_repeats
print('Running Python analyze_mtDNA_repeats on half-genome...')
analyze_mtDNA_repeats(str(fasta), pos, ref, alt, output_path=str(tmp), output_prefix=outpref, debug=True)
print('Python finished.')

# compare outputs
rA = tmp / f"01KP.{pos}.{ref}.txt"
pyA = tmp / f"pos{pos}_{ref}to{alt}" / f"01KP.{pos}.{ref}.txt"

print('\nFiles:')
print('R:', rA.exists(), rA)
print('PY:', pyA.exists(), pyA)

def load_tab(p):
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p, sep='\t', dtype=str)

r_df = load_tab(rA)
py_df = load_tab(pyA)
print('Rows R:', len(r_df), 'Rows PY:', len(py_df))

if not r_df.empty or not py_df.empty:
    all_cols = sorted(set(r_df.columns).union(set(py_df.columns)))
    for c in all_cols:
        if c not in r_df.columns: r_df[c] = ""
        if c not in py_df.columns: py_df[c] = ""
    r_s = r_df[all_cols].astype(str).sort_values(by=all_cols).reset_index(drop=True)
    py_s = py_df[all_cols].astype(str).sort_values(by=all_cols).reset_index(drop=True)
    try:
        pd.testing.assert_frame_equal(r_s, py_s, check_dtype=False)
        print('DataFrames match exactly')
    except AssertionError as e:
        print('DataFrames differ')
        # print first mismatch
        max_rows = max(len(r_s), len(py_s))
        if len(r_s) < max_rows:
            r_s = pd.concat([r_s, pd.DataFrame([[''] * len(all_cols)] * (max_rows - len(r_s)), columns=all_cols)], ignore_index=True)
        if len(py_s) < max_rows:
            py_s = pd.concat([py_s, pd.DataFrame([[''] * len(all_cols)] * (max_rows - len(py_s)), columns=all_cols)], ignore_index=True)
        mismatch = (r_s.values != py_s.values)
        idx = int(mismatch.argmax())
        row = idx // mismatch.shape[1]
        col = idx % mismatch.shape[1]
        print('First mismatch at row', row, 'column', all_cols[col])
        print('R:', r_s.iloc[row, col])
        print('PY:', py_s.iloc[row, col])

print('\nDone')
