# REPEATABILITY — Python port

This directory contains a Python reimplementation of the R functions from
`scripts/01.RefAltRepeatability.Function.Rmd`.

Structure:
- `src/io.py` — reading FASTA and CSV utilities
- `src/analysis.py` — core motif and repeat-finding logic
- `src/viz.py` — plotting helpers (minimal)
- `src/main.py` — example runner that replicates the R workflow for pos=8473
- `tests/test_analysis.py` — minimal unit test

Requirements: see `requirements.txt`.

Run the example (from this folder):

```powershell
# create venv, install requirements, then
python -m src.main
```

Optimizations suggested & implemented:
- vectorized k-mer comparison using NumPy (implemented)
- optional numba JIT for tight loops (dependency added to `requirements.txt`)
- future work: multiprocessing / Dask for broad parallelism, suffix-array indexing for many queries, Cython/C extension for ultimate speed

How to enable optimizations:
- In analysis functions you can pass `use_numpy=True` to enable NumPy vectorized k-mer/Hamming search.
- If `numba` is installed, certain hot functions can be JIT-compiled (automatic where available).

