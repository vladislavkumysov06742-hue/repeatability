import pandas as pd
from pathlib import Path

tmp = Path(__file__).resolve().parents[1] / '.pytest_tmp' / 'debug_run_tmp'
pos = 100
ref = 'A'
alt = 'G'
py_folder = tmp / f'pos{pos}_{ref}to{alt}'

files = [
    (tmp / f'01KP.{pos}.{ref}.txt', py_folder / f'01KP.{pos}.{ref}.txt'),
    (tmp / f'01KP.{pos}.{alt}.txt', py_folder / f'01KP.{pos}.{alt}.txt')
]

for r_file, py_file in files:
    print('\nComparing', r_file, 'vs', py_file)
    if not r_file.exists():
        print('Missing R file', r_file)
        continue
    if not py_file.exists():
        print('Missing PY file', py_file)
        continue
    r = pd.read_csv(r_file, sep='\t', dtype=str)
    py = pd.read_csv(py_file, sep='\t', dtype=str)
    all_cols = sorted(set(r.columns).union(py.columns))
    for c in all_cols:
        if c not in r.columns: r[c] = ''
        if c not in py.columns: py[c] = ''
    r = r[all_cols].astype(str)
    py = py[all_cols].astype(str)
    set_r = set(tuple(x) for x in r.values.tolist())
    set_py = set(tuple(x) for x in py.values.tolist())
    only_r = set_r - set_py
    only_py = set_py - set_r
    print('rows R:', len(r), 'rows PY:', len(py))
    print('only in R:', len(only_r), 'only in PY:', len(only_py))
    if len(only_r) > 0:
        print('\nSample rows only in R:')
        for i, row in enumerate(list(only_r)[:10]):
            print({col: val for col, val in zip(all_cols, row)})
            # check whether these R-only rows were present in Python raw_all_repeats (before filtering)
            dbg_folder = tmp / 'debug' / f'pos{pos}_{ref}to{alt}'
            raw_file = dbg_folder / 'my_test_repeat_raw_all_repeats.csv'
            if raw_file.exists():
                raw = pd.read_csv(raw_file, dtype=str)
                raw_cols = sorted(set(raw.columns))
                for i, row in enumerate(list(only_r)[:10]):
                    row_dict = {col: val for col, val in zip(all_cols, row)}
                    # build boolean mask by matching repeat.start and repeat.seq if available
                    mask = pd.Series([True] * len(raw))
                    if 'repeat.start' in raw.columns and row_dict.get('repeat.start') is not None:
                        mask = mask & (raw['repeat.start'].astype(str) == str(row_dict.get('repeat.start')))
                    if 'repeat.seq' in raw.columns and row_dict.get('repeat.seq') is not None:
                        mask = mask & (raw['repeat.seq'].astype(str) == str(row_dict.get('repeat.seq')))
                    found = raw[mask]
                    print(' R-only row present in Python raw_all_repeats?', 'YES' if not found.empty else 'NO')
    if len(only_py) > 0:
        print('\nSample rows only in PY:')
        for i, row in enumerate(list(only_py)[:10]):
            print({col: val for col, val in zip(all_cols, row)})

# Also compare Python debug non_nested vs final py outputs for sanity
dbg_folder = tmp / 'debug' / f'pos{pos}_{ref}to{alt}'
if dbg_folder.exists():
    print('\nPython debug files in', dbg_folder)
    for name in ['my_test_repeat_raw_all_repeats.csv', 'my_test_repeat_filtered_self_overlaps.csv', 'my_test_repeat_non_nested.csv']:
        p = dbg_folder / name
        print(name, 'exists' if p.exists() else 'MISSING')
else:
    print('\nNo Python debug folder found')

print('\nDone')
