from pathlib import Path
import pandas as pd

BASE = Path('.').resolve() / '.pytest_tmp'
if not BASE.exists():
    print('No .pytest_tmp found in', Path('.').resolve())
    raise SystemExit(1)

# Find the latest test_r_and_python_equivalence* directory
cands = [p for p in BASE.iterdir() if p.is_dir() and p.name.startswith('test_r_and_python_equivalence')]
if not cands:
    print('No test_r_and_python_equivalence* dirs found in', BASE)
    raise SystemExit(1)

latest = max(cands, key=lambda p: p.stat().st_mtime)
print('Using tmp dir:', latest)

r_file = latest / '01KP.100.A.txt'
py_file = latest / 'pos100_AtoG' / '01KP.100.A.txt'
print('R file:', r_file)
print('Py file:', py_file)

if not r_file.exists():
    print('R file missing'); raise SystemExit(2)
if not py_file.exists():
    print('Py file missing'); raise SystemExit(3)

r_df = pd.read_csv(r_file, sep='\t', dtype=str)
py_df = pd.read_csv(py_file, sep='\t', dtype=str)

# Ensure same length for row-wise comparison
min_rows = min(len(r_df), len(py_df))
print(f'Rows: R={len(r_df)}, Py={len(py_df)}, comparing first {min_rows} rows')

out = []
for i in range(min_rows):
    r_m = r_df.at[i, 'motif.seq'] if 'motif.seq' in r_df.columns else ''
    p_m = py_df.at[i, 'motif.seq'] if 'motif.seq' in py_df.columns else ''
    if r_m != p_m:
        # find first differing index
        diff_pos = None
        maxlen = max(len(r_m), len(p_m))
        for k in range(maxlen):
            a = r_m[k] if k < len(r_m) else None
            b = p_m[k] if k < len(p_m) else None
            if a != b:
                diff_pos = k
                break
        row = {
            'index': i,
            'repeat.start_R': r_df.at[i, 'repeat.start'] if 'repeat.start' in r_df.columns else '',
            'repeat.end_R': r_df.at[i, 'repeat.end'] if 'repeat.end' in r_df.columns else '',
            'motif_R': r_m,
            'motif_Py': p_m,
            'first_diff_pos': diff_pos,
        }
        out.append(row)
        if len(out) >= 20:
            break

if not out:
    print('No differences found in compared rows')
else:
    import json
    print(json.dumps(out, indent=2))
