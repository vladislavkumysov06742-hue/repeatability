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
    # Load full reference mtDNA sequence
    repo_root = Path(__file__).resolve().parents[2]
    ref_genome_file = repo_root / "data" / "1_raw" / "Homo_sapients.mtDNA.fasta"
    assert ref_genome_file.exists(), f"Reference genome file not found: {ref_genome_file}"
    
    # Read full reference sequence
    full_seq = ""
    with open(ref_genome_file) as f:
        for line in f:
            if not line.startswith(">"):
                full_seq += line.strip()
    
    # Create test fasta file with full sequence
    fasta = tmp_path / "test_full_mtDNA.fasta"
    fasta.write_text(">human_mtDNA\n" + full_seq + "\n")
    
    # Path to R functions file in repository
    r_functions = repo_root / "scripts" / "02.RefAltRepeatability.Functions.R"
    assert r_functions.exists(), f"R functions file not found: {r_functions}"
    
    # Test the three biologically real positions in range 8000-9000
    # These are the positions where we have actual mtDNA mutations
    test_positions = [
        (8251, "G", "A"),  # pos8251_GtoA
        (8472, "C", "T"),  # pos8472_CtoT
        (8473, "T", "C"),  # pos8473_TtoC
    ]
    
    for pos, ref, alt in test_positions:
        # Verify position is within genome range
        assert 1 <= pos <= len(full_seq), f"Position {pos} is outside genome range (1-{len(full_seq)})"
        
        # Verify the reference nucleotide matches the sequence at this position
        # Convert 1-based position to 0-based index
        actual_nuc = full_seq[pos - 1]
        assert actual_nuc == ref, f"Position {pos}: expected {ref}, found {actual_nuc}"
        
        outpref = "my_mtDNA_repeat"
        
        # Create R runner script for this position
        r_runner = tmp_path / f"run_r_analysis_pos{pos}.R"
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
            assert r_file.exists(), f"R output file missing for pos {pos}, allele {allele}: {r_file}"
            assert py_file.exists(), f"Python output file missing for pos {pos}, allele {allele}: {py_file}"
            
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
