from Bio import SeqIO
import pandas as pd


def read_fasta_sequence(path):
    """Read the first sequence from a FASTA file and return as plain string."""
    rec = SeqIO.read(path, "fasta")
    return str(rec.seq)


def read_csv(path, **kwargs):
    return pd.read_csv(path, **kwargs)


def write_csv(df, path, **kwargs):
    df.to_csv(path, index=False, **kwargs)
