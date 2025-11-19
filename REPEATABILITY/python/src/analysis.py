"""Core analysis: motif extraction, approximate repeat finding, filtering and repeatability summary."""
from typing import List, Dict
import math
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import os
import numpy as np
from numpy.lib.stride_tricks import as_strided
import logging
import csv

try:
    import numba
    NUMBA_AVAILABLE = True
except Exception:
    NUMBA_AVAILABLE = False

logger = logging.getLogger(__name__)


def get_motif_around_position(seq: str, pos: int, left_flank: int, right_flank: int) -> str:
    # pos is 1-based to match R; convert to 0-based
    start = max(pos - left_flank - 1, 0)
    end = min(pos + right_flank - 1, len(seq) - 1)
    return seq[start:end + 1]


def hamming_distance(a: str, b: str) -> int:
    # assume equal length
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def _seq_to_int_array(seq: str) -> np.ndarray:
    """Map sequence string to small int array (A:0,C:1,G:2,T:3,others:4)."""
    mp = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3, ord('a'):0, ord('c'):1, ord('g'):2, ord('t'):3}
    arr = np.frombuffer(seq.encode('ascii', 'replace'), dtype=np.uint8)
    # map bytes to our small ints
    out = np.full(arr.shape, 4, dtype=np.uint8)
    for k, v in mp.items():
        out[arr == k] = v
    return out


def _sliding_window_view(arr: np.ndarray, k: int) -> np.ndarray:
    """Return a 2D view (n-k+1, k) over 1D array using strides."""
    n = arr.shape[0]
    if k > n:
        return np.empty((0, k), dtype=arr.dtype)
    shape = (n - k + 1, k)
    strides = (arr.strides[0], arr.strides[0])
    return as_strided(arr, shape=shape, strides=strides)


def find_approximate_repeats_np(seq: str, motif: str, max_mismatch: int) -> pd.DataFrame:
    """NumPy-vectorized implementation of approximate repeats using sliding window and vector ops.

    Returns same columns as the pure-Python implementation.
    """
    k = len(motif)
    n = len(seq)
    if k > n:
        return pd.DataFrame(columns=["repeat.seq", "repeat.start", "repeat.end", "repeat.hamming.distance"])
    seq_arr = _seq_to_int_array(seq)
    motif_arr = _seq_to_int_array(motif)
    windows = _sliding_window_view(seq_arr, k)
    # compute mismatches per window
    # boolean matrix where True indicates mismatch
    mismatches = windows != motif_arr
    distances = mismatches.sum(axis=1)
    matches_idx = np.nonzero(distances <= max_mismatch)[0]
    rows = []
    if matches_idx.size > 0:
        # reconstruct strings for matched kmers
        # to avoid creating all kmers, build from original string slicing
        for i in matches_idx.tolist():
            rows.append({
                "repeat.seq": seq[i:i + k],
                "repeat.start": i + 1,
                "repeat.end": i + k,
                "repeat.hamming.distance": int(distances[i])
            })
    return pd.DataFrame(rows)


def find_approximate_repeats(seq: str, motif: str, max_mismatch: int, use_numpy: bool = False, use_numba: bool = False) -> pd.DataFrame:
    """Find all k-mers in seq that differ from motif by <= max_mismatch.

    Parameters
    - use_numpy: if True, use a NumPy-vectorized implementation (much faster for long sequences).
    - use_numba: reserved for numba-based implementation (if available); currently unused when numba is not installed.
    """
    if use_numpy:
        try:
            return find_approximate_repeats_np(seq, motif, max_mismatch)
        except Exception as e:
            logger.warning("NumPy accelerated search failed, falling back to Python implementation: %s", e)
    # fallback pure-Python
    k = len(motif)
    n = len(seq)
    if k > n:
        return pd.DataFrame(columns=["repeat.seq", "repeat.start", "repeat.end", "repeat.hamming.distance"])
    kmers = [seq[i:i + k] for i in range(0, n - k + 1)]
    distances = [hamming_distance(motif, s) for s in kmers]
    matches = [i for i, d in enumerate(distances) if d <= max_mismatch]
    rows = []
    for i in matches:
        rows.append({
            "repeat.seq": kmers[i],
            "repeat.start": i + 1,  # 1-based
            "repeat.end": i + k,
            "repeat.hamming.distance": distances[i]
        })
    return pd.DataFrame(rows)


def _process_task_worker(params):
    """Top-level worker function for ProcessPoolExecutor. Receives a tuple of primitive args.

    params: (allele, seq_allele, motif_length, left_flank, right_flank, pos, ref_nuc, use_numpy, use_numba)
    Returns a pandas DataFrame (possibly empty).
    """
    allele, seq_allele, motif_length, left_flank, right_flank, pos, ref_nuc, use_numpy, use_numba = params

    # build motif around position using provided flanks
    motif_seq = get_motif_around_position(seq_allele, pos, left_flank, right_flank)
    motif_string = motif_seq
    this_length = len(motif_string)
    max_mismatch = math.floor(0.2 * this_length)

    # find repeats (numpy or python implementation)
    repeats_df = find_approximate_repeats(seq_allele, motif_seq, max_mismatch,
                                          use_numpy=use_numpy, use_numba=use_numba)

    if repeats_df.empty:
        return pd.DataFrame()

    # annotate and normalize columns
    repeats_df = repeats_df.copy()
    repeats_df["motif.seq"] = motif_string
    repeats_df["motif.length"] = this_length
    repeats_df["motif.start"] = max(pos - left_flank, 1)
    repeats_df["motif.end"] = min(pos + right_flank, len(seq_allele))
    repeats_df["nuc"] = allele
    repeats_df["pos"] = pos
    repeats_df["RefAlt"] = "Ref" if allele == ref_nuc else "Alt"

    cols = ["pos", "nuc", "RefAlt", "motif.seq", "motif.length", "motif.start", "motif.end",
            "repeat.seq", "repeat.start", "repeat.end", "repeat.hamming.distance"]
    available = [c for c in cols if c in repeats_df.columns]
    return repeats_df[available]


def get_repeatability_table(seq: str, pos: int, ref_nuc: str, alt_nuc: str,
                            max_flank_left: int = 20, max_flank_right: int = 20,
                            min_length: int = 5, max_length: int = 41,
                            use_numpy: bool = False,
                            use_numba: bool = False,
                            use_parallel: bool = False,
                            n_workers: int | None = None) -> pd.DataFrame:
    """Compute repeats for reference and alt alleles similar to the R routine.

    New options:
    - use_numpy: use NumPy-vectorized k-mer/hamming implementation
    - use_numba: use numba JIT (if available) for inner loops
    - use_parallel: parallelize motif searches using ProcessPoolExecutor
    """
    results = []

    def _allele_sequence(allele_char: str) -> str:
        if allele_char == ref_nuc:
            return seq
        seq_list = list(seq)
        seq_list[pos - 1] = allele_char
        return "".join(seq_list)

    # Build list of tasks so we can parallelize if requested
    tasks = []
    for allele in (ref_nuc, alt_nuc):
        seq_allele = _allele_sequence(allele)
        for motif_length in range(min_length, max_length + 1):
            for left_flank in range(0, motif_length):
                right_flank = motif_length - left_flank - 1
                if left_flank <= max_flank_left and right_flank <= max_flank_right:
                    # include parameters needed by top-level worker
                    tasks.append((allele, seq_allele, motif_length, left_flank, right_flank, pos, ref_nuc, use_numpy, (use_numba and NUMBA_AVAILABLE)))

    if use_parallel:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        results_dfs = []
        with ProcessPoolExecutor(max_workers=n_workers) as exe:
            futures = {exe.submit(_process_task_worker, task): task for task in tasks}
            for fut in as_completed(futures):
                try:
                    df = fut.result()
                    if not df.empty:
                        results_dfs.append(df)
                except Exception as e:
                    logger.exception("Task failed: %s", e)
        results = results_dfs
    else:
        for task in tasks:
            df = _process_task_worker(task)
            if not df.empty:
                results.append(df)

    if len(results) == 0:
        return pd.DataFrame(columns=["pos", "nuc", "RefAlt", "motif.seq", "motif.length", "motif.start", "motif.end",
                                     "repeat.seq", "repeat.start", "repeat.end", "repeat.hamming.distance"])
    return pd.concat(results, ignore_index=True)


def analyze_mtDNA_repeats(fasta_path: str,
                          pos: int,
                          ref_nuc: str,
                          alt_nuc: str,
                          max_flank_left: int = 20,
                          max_flank_right: int = 20,
                          min_length: int = 5,
                          mismatch_percent: float = 0.2,
                          use_numpy: bool = False,
                          use_numba: bool = False,
                          use_parallel: bool = False,
                          n_workers: int | None = None,
                          major_arc_start: int = 5798,
                          major_arc_end: int = 16568,
                          output_path: str = None,
                          output_prefix: str = "my_mtDNA_repeat") -> Dict[str, Path]:
    """High-level pipeline mirroring the R `analyze_mtDNA_repeats` function.

    Reads fasta, finds repeats for Ref and Alt, filters self-overlaps and nested repeats,
    computes EffectiveLength, writes two CSVs into output folder and returns paths.
    """
    # Read sequence
    rec = SeqIO.read(str(fasta_path), "fasta")
    seq = str(rec.seq)

    # Ensure that the nucleotide at `pos` equals the provided reference allele (match R behavior)
    if pos < 1 or pos > len(seq):
        raise ValueError(f"pos {pos} out of range for sequence length {len(seq)}")
    if seq[pos - 1] != ref_nuc:
        seq_list = list(seq)
        seq_list[pos - 1] = ref_nuc
        seq = "".join(seq_list)

    # Compute raw repeats table
    table = get_repeatability_table(seq, pos, ref_nuc, alt_nuc,
                                     max_flank_left=max_flank_left, max_flank_right=max_flank_right,
                                     min_length=min_length, max_length=(max_flank_left + max_flank_right + 1),
                                     use_numpy=use_numpy, use_numba=use_numba, use_parallel=use_parallel, n_workers=n_workers)

    if table.empty:
        # prepare output folder even if empty
        folder = Path(output_path) if output_path else Path(".")
        folder = folder / f"pos{pos}_{ref_nuc}to{alt_nuc}"
        folder.mkdir(parents=True, exist_ok=True)
        all_out = folder / f"{output_prefix}_all_repeats.csv"
        summary_out = folder / f"{output_prefix}_major_arc_summary_top5.csv"
        pd.DataFrame().to_csv(all_out, index=False)
        pd.DataFrame().to_csv(summary_out, index=False)
        return {"all_repeats": all_out, "summary_top5_major_arc": summary_out, "output_folder": folder}

    # Filter self overlaps and compute EffectiveLength
    filtered = filter_self_overlaps(table)
    filtered = compute_effective_length(filtered)

    # Remove nested repeats
    non_nested = remove_nested_repeats(filtered)

    # Prepare major-arc summary (top5 effective lengths) similar to R script
    maj = non_nested[(non_nested["repeat.start"] >= major_arc_start) & (non_nested["repeat.end"] <= major_arc_end)].copy()
    if maj.empty:
        summary_df = pd.DataFrame()
    else:
        top5 = maj.groupby(["RefAlt", "EffectiveLength"]).size().reset_index(name="Count")
        top5 = top5[top5["EffectiveLength"] > 0]
        # For each RefAlt take top 5 EffectiveLength by descending
        top5 = top5.sort_values(["RefAlt", "EffectiveLength"], ascending=[True, False])
        top5 = top5.groupby("RefAlt").head(5)
        # Build wide summary similar to R
        efflen_top = sorted(top5["EffectiveLength"].unique(), reverse=True)[:5]
        summary_rows = []
        for eff in efflen_top:
            ref_count = int(((top5["RefAlt"] == "Ref") & (top5["EffectiveLength"] == eff)).sum())
            alt_count = int(((top5["RefAlt"] == "Alt") & (top5["EffectiveLength"] == eff)).sum())
            summary_rows.append({"EffectiveLength": int(eff), "Ref_Count": ref_count, "Alt_Count": alt_count})
        summary_df = pd.DataFrame(summary_rows)

    # Create output folder
    folder = Path(output_path) if output_path else Path(".")
    folder = folder / f"pos{pos}_{ref_nuc}to{alt_nuc}"
    folder.mkdir(parents=True, exist_ok=True)
    all_out = folder / f"{output_prefix}_all_repeats.csv"
    summary_out = folder / f"{output_prefix}_major_arc_summary_top5.csv"

    # Write outputs
    non_nested.to_csv(all_out, index=False)
    summary_df.to_csv(summary_out, index=False)

    # Also write Greedy-style per-allele tab-delimited outputs: 01KP.<pos>.<nuc>.txt
    if not non_nested.empty:
        # add lowercase effective.length column for compatibility
        non_nested = non_nested.copy()
        if "EffectiveLength" in non_nested.columns:
            non_nested["effective.length"] = non_nested["EffectiveLength"]
        else:
            non_nested["effective.length"] = non_nested["motif.length"] - non_nested["repeat.hamming.distance"]
        # ensure effective.length is numeric and sort deterministically to match R's output ordering
        non_nested = non_nested.copy()
        # convert to numeric if possible
        try:
            non_nested["effective.length"] = pd.to_numeric(non_nested["effective.length"])
        except Exception:
            # leave as-is if conversion fails
            pass
        # R sorts per-allele output by effective.length descending only; replicate that to match R's ordering
        if "effective.length" in non_nested.columns:
            non_nested = non_nested.sort_values(by=["effective.length"], ascending=[False])
        for allele in non_nested["nuc"].unique():
            sub = non_nested[non_nested["nuc"] == allele]
            if sub.empty:
                continue
            cols = ["pos", "nuc", "motif.seq", "motif.length", "motif.start", "motif.end",
                    "repeat.seq", "repeat.start", "repeat.end", "repeat.hamming.distance", "effective.length"]
            out_file = folder / f"01KP.{pos}.{allele}.txt"
            # Build an output DataFrame that always contains the R-expected columns in order.
            # Fill missing columns with empty strings so files are produced deterministically.
            out_df = pd.DataFrame()
            for c in cols:
                if c in sub.columns:
                    out_df[c] = sub[c].astype(str)
                else:
                    out_df[c] = [""] * len(sub)
            out_df.to_csv(out_file, sep="\t", index=False)

    return {"all_repeats": all_out, "summary_top5_major_arc": summary_out, "output_folder": folder}


def filter_self_overlaps(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    return df[(df["motif.end"] < df["repeat.start"]) | (df["repeat.end"] < df["motif.start"])].copy()


def compute_effective_length(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df["EffectiveLength"] = df["motif.length"] - df["repeat.hamming.distance"]
    return df


def remove_nested_repeats(df: pd.DataFrame) -> pd.DataFrame:
    """Greedy removal of repeats nested inside previously kept longer repeats per allele."""
    if df.empty:
        return df
    kept = []
    for allele in df["nuc"].unique():
        sub = df[df["nuc"] == allele].copy()
        sub = sub.sort_values(by="EffectiveLength", ascending=False)
        kept_intervals = []  # list of (start,end)
        for _, row in sub.iterrows():
            cs, ce = int(row["repeat.start"]), int(row["repeat.end"])
            # check if inside any kept interval
            is_nested = any(cs >= ks and ce <= ke for ks, ke in kept_intervals)
            if not is_nested:
                kept.append(row)
                kept_intervals.append((cs, ce))
    if len(kept) == 0:
        return pd.DataFrame(columns=df.columns)
    return pd.DataFrame(kept).reset_index(drop=True)


def summarize_major_arc(df: pd.DataFrame, major_arc_start: int = 5798, major_arc_end: int = 16568) -> pd.DataFrame:
    if df.empty:
        return df
    maj = df[(df["motif.start"] >= major_arc_start) & (df["motif.end"] <= major_arc_end) &
             (df["repeat.start"] >= major_arc_start) & (df["repeat.end"] <= major_arc_end)].copy()
    if maj.empty:
        return maj
    max_eff = int(maj["EffectiveLength"].max())
    rows = []
    for eff in range(max_eff, max(max_eff - 5, 0) - 1, -1):
        ref_count = int(((maj["EffectiveLength"] == eff) & (maj["nuc"] == "T")).sum())
        alt_count = int(((maj["EffectiveLength"] == eff) & (maj["nuc"] == "C")).sum())
        rows.append({"EffectiveLength": eff, "Ref_Count": ref_count, "Alt_Count": alt_count})
    return pd.DataFrame(rows)
