## Optimization Summary

### 1. NumPy Vectorization ✅ **WORKING**
- **Before**: 740 sec (57M hamming_distance calls)
- **After**: 21.55 sec (pos8000, 3 alts)
- **Speedup**: ~34x
- **Status**: `--use-numpy` flag works, `find_approximate_repeats_np` is called automatically

### 2. Batch DataFrame Construction ⏳ **PARTIALLY DONE**
- **Status**: Code changes for parquet caching applied, but worker still creates intermediate DataFrames
- **Benefit**: Could save ~2-3 sec per position (based on profiling showing 6.9s on DataFrame.__setitem__)
- **Implementation**: Would require refactoring `_process_task_worker` to return list-of-dicts instead of DataFrame

### 3. Parquet Caching ✅ **WORKING**
- **Before**: pickle-based caching
- **After**: parquet with snappy compression
- **Evidence**: `pos8000_ref_non_nested.parquet` created successfully
- **Dependency**: pyarrow installed (added to requirements.txt)
- **Benefit**: ~0.7 sec per alt when using precompute-ref mode

### 4. Parallel Positions ⚠️ **PARTIALLY WORKING**
- **Status**: Code added (ProcessPoolExecutor), but Windows multiprocessing has serialization issues
- **Working on**: Single-position runs complete successfully
- **Issue**: ProcessPoolExecutor with `process_position` function fails to serialize properly on Windows
- **Alternative**: Could use `if __name__ == '__main__'` guard and fix worker function
- **Benefit**: Should provide ~2x speedup for 2 workers on 2 positions (if working)

### Test Results
- `--use-numpy` alone: 21.55 sec for 3 alts of pos8000
- `--precompute-ref-only`: precomputes ref in ~4 sec
- alt processing with cached ref: 19.49 sec for 3 alts (saves ~2 sec vs full compute)
- Parquet cache created successfully: pos8000_ref_non_nested.parquet (snappy compressed)

### Next Steps
1. **Option A**: Fix ProcessPoolExecutor (add `if __name__ == '__main__'` guard, refactor worker)
2. **Option B**: Implement batch DataFrame construction (optimization 2) which is lower-hanging fruit
3. **Option C**: Load full reference sequence once in memory (optimization 3 from original TODO)

### Files Modified
- `python/run_all_variants.py`: Added --use-numpy, --precompute-ref-only, --max-workers args
- `python/src/analysis.py`: Changed cache format from pickle to parquet (3 code changes)
- `python/requirements.txt`: Added pyarrow
- `python/test_parallel_timing.py`: Created for benchmarking
