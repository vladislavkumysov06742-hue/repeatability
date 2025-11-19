import subprocess
import shutil
from pathlib import Path
import sys
import pandas as pd
import pytest

from src.analysis import analyze_mtDNA_repeats


def is_rscript_available():
    return shutil.which("Rscript") is not None


@pytest.mark.skipif(not is_rscript_available(), reason="Rscript not available on PATH")
def test_r_and_python_equivalence(tmp_path):
    # Prepare synthetic fasta
    seq = ("GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCACTTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAATCTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATACCCCGAACCAACCAAACCCCAAAGACACCCCCCACAGTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACAAATAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGTTCCAGTGAGTTCACCCTCTAAATCACCACGATCAAAAGGAACAAGCATCAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCCACACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGATACCCCACTATGCTTAGCCCTAAACCTCAACAGTTAAATCAACAAAACTGCTCGCCAGAACACTACGAGCCACAGCTTAAAACTCAAAGGACCTGGCGGTGCTTCATATCCCTCTAGAGGAGCCTGTTCTGTAATCGATAAACCCCGATCAACCTCACCACCTCTTGCTCAGCCTATATACCGCCATCTTCAGCAAACCCTGATGAAGGCTACAAAGTAAGCGCAAGTACCCACGTAAAGACGTTAGGTCAAGGTGTAGCCCATGAGGTGGCAAGAAATGGGCTACATTTTCTACCCCAGAAAACTACGATAGCCCTTATGAAACTTAAGGGTCGAAGGTGGATTTAGCAGTAAACTAAGAGTAGAGTGCTTAGTTGAACAGGGCCCTGAAGCGCGTACACACCGCCCGTCACCCTCCTCAAGTATACTTCAAAGGACATTTAACTAAAACCCCTACGCATTTATATAGAGGAGACAAGTCGTAACATGGTAAGTGTACTGGAAAGTGCACTTGGACGAACCAGAGTGTAGCTTAACACAAAGCACCCAACTTACACTTAGGAGATTTCAACTTAACTTGACCGCTCTGAGCTAAACCTAGCCCCAAACCCACTCCACCTTACTACCAGACAACCTTAGCCAAACCATTTACCCAAATAAAGTATAGGCGATAGAAATTGAAACCTGGCGCAATAGATATAGTACCGCAAGGGAAAGATGAAAAATTATAACCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAATTAACTAGAAATAACTTTGCAAGGAGAGCCAAAGCTAAGACCCCCGAAACCAGACGAGCTACCTAAGAACAGCTAAAAGAGCACACCCGTCTATGTAGCAAAATAGTGGGAAGATTTATAGGTAGAGGCGACAAACCTACCGAGCCTGGTGATAGCTGGTTGTCCAAGAT")[:2000]
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">test\n" + seq + "\n")

    pos = 100
    ref = "A"
    alt = "G"
    outpref = "my_test_repeat"

    # Path to R functions file in repository
    repo_root = Path(__file__).resolve().parents[2]
    r_functions = repo_root / "scripts" / "02.RefAltRepeatability.Functions.R"
    assert r_functions.exists(), f"R functions file not found: {r_functions}"

    # Create small R runner script in tmp_path
    r_runner = tmp_path / "run_r_analysis.R"
    r_runner.write_text(f"source('{r_functions.as_posix()}')\n" \
                      f"analyze_mtDNA_repeats(fasta_path = '{fasta.as_posix()}', pos = {pos}, ref_nuc = '{ref}', alt_nuc = '{alt}', output_prefix = '{outpref}')\n")

    # Run R
    subprocess.run(["Rscript", r_runner.as_posix()], check=True, cwd=tmp_path)

    # Run Python, writing outputs into tmp_path as output_path
    analyze_mtDNA_repeats(str(fasta), pos, ref, alt, output_path=str(tmp_path), output_prefix=outpref)

    # Compare per-allele 01KP files
    for allele in (ref, alt):
        r_file = tmp_path / f"01KP.{pos}.{allele}.txt"
        py_file = tmp_path / f"pos{pos}_{ref}to{alt}" / f"01KP.{pos}.{allele}.txt"
        assert r_file.exists(), f"R output file missing: {r_file}"
        assert py_file.exists(), f"Python output file missing: {py_file}"

        # Read as DataFrame (tab-delimited)
        r_df = pd.read_csv(r_file, sep="\t", dtype=str)
        py_df = pd.read_csv(py_file, sep="\t", dtype=str)

        # Normalize column names
        r_df.columns = [c.strip() for c in r_df.columns]
        py_df.columns = [c.strip() for c in py_df.columns]

        # To avoid false failures due to different row ordering, compare as sets:
        # align columns, fill missing cols with empty string, then sort by all columns.
        all_cols = sorted(set(r_df.columns).union(set(py_df.columns)))
        for c in all_cols:
            if c not in r_df.columns:
                r_df[c] = ""
            if c not in py_df.columns:
                py_df[c] = ""

        r_df_sorted = r_df[all_cols].astype(str).sort_values(by=all_cols).reset_index(drop=True)
        py_df_sorted = py_df[all_cols].astype(str).sort_values(by=all_cols).reset_index(drop=True)

        # Compare contents ignoring dtypes and original index
        pd.testing.assert_frame_equal(r_df_sorted, py_df_sorted, check_like=True, check_dtype=False)
