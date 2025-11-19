"""Example runner that replicates parts of the R workflow for position 8473."""
from pathlib import Path
from .io import read_fasta_sequence, write_csv
from .analysis import (
    get_repeatability_table, filter_self_overlaps, compute_effective_length,
    remove_nested_repeats, summarize_major_arc
)
import pandas as pd


def run_example():
    base = Path(__file__).resolve().parents[1]
    fasta = base / ".." / "data" / "1_raw" / "Homo_sapients.mtDNA.fasta"
    fasta = fasta.resolve()
    seq = read_fasta_sequence(str(fasta))
    pos = 8473
    ref = "T"
    alt = "C"
    print("Computing repeats (this may take a few seconds)...")
    table = get_repeatability_table(seq, pos, ref, alt)
    table = filter_self_overlaps(table)
    table = compute_effective_length(table)
    non_nested = remove_nested_repeats(table)
    summary = summarize_major_arc(non_nested)
    out_dir = base / "output"
    out_dir.mkdir(exist_ok=True)
    write_csv(non_nested, str(out_dir / "non_nested_repeats.csv"))
    write_csv(summary, str(out_dir / "summary_repeatability.csv"))
    print("Wrote results to:", out_dir)


if __name__ == "__main__":
    run_example()
