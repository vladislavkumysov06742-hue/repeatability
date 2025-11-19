import pandas as pd
import numpy as np
from src.analysis import (
    get_motif_around_position,
    find_approximate_repeats,
    get_repeatability_table,
    remove_nested_repeats,
    compute_effective_length,
    summarize_major_arc,
)


def test_get_motif_exact():
    seq = "ACGTACGT"
    # pos=3 (G), left=1, right=1 -> positions 2..4 => "CGT"
    motif = get_motif_around_position(seq, pos=3, left_flank=1, right_flank=1)
    assert motif == "CGT"


def test_find_approximate_repeats_counts():
    seq = "AAAAACCCCCAAAAACCCCC"
    motif = "AAAAA"
    df = find_approximate_repeats(seq, motif, max_mismatch=0)
    # Expect two perfect matches (positions 1 and 11)
    assert list(df['repeat.start']) == [1, 11]


def test_get_repeatability_table_columns_and_alleles():
    seq = "AAAAACCCCCAAAAACCCCC"
    pos = 3
    ref = "A"
    alt = "G"
    table = get_repeatability_table(seq, pos, ref, alt, max_flank_left=2, max_flank_right=2, min_length=3, max_length=3)
    # Should return a DataFrame with required columns and contain both alleles
    assert not table.empty
    for col in ["nuc", "motif.seq", "motif.length", "repeat.seq", "repeat.start"]:
        assert col in table.columns
    assert set(table['nuc'].unique()).issuperset({ref, alt})


def test_remove_nested_repeats_behavior():
    # Construct synthetic repeats: one long and one nested inside it for same allele
    df = pd.DataFrame({
        'nuc': ['A', 'A', 'C'],
        'repeat.start': [1, 3, 20],
        'repeat.end': [10, 5, 25],
        'motif.length': [10, 3, 6],
        'repeat.hamming.distance': [0, 0, 0],
    })
    df = compute_effective_length(df)
    kept = remove_nested_repeats(df)
    # For allele 'A' only the longer interval (1-10) should remain
    a_kept = kept[kept['nuc'] == 'A']
    assert len(a_kept) == 1
    assert a_kept.iloc[0]['repeat.start'] == 1


def test_summarize_major_arc_counts():
    # Create DF with EffectiveLength and repeat positions inside major arc
    df = pd.DataFrame({
        'nuc': ['T', 'C', 'T'],
        'motif.start': [5800, 6000, 7000],
        'motif.end': [5804, 6004, 7004],
        'repeat.start': [5800, 6000, 7000],
        'repeat.end': [5804, 6004, 7004],
        'EffectiveLength': [16, 15, 14]
    })
    summary = summarize_major_arc(df, major_arc_start=5798, major_arc_end=16568)
    # summary should have rows for EffectiveLength values (descending)
    assert 'Ref_Count' in summary.columns or 'Alt_Count' in summary.columns or 'Ref_Count' in summary.columns
    assert not summary.empty
