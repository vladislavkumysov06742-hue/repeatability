import sys
from pathlib import Path

try:
    import pandas as pd
except Exception as e:
    print('Pandas import error:', e)
    sys.exit(2)

try:
    from src.analysis import analyze_mtDNA_repeats
except Exception as e:
    print('Import error from src.analysis:', e)
    sys.exit(2)

base = Path(__file__).resolve().parents[0]
# input fasta (relative to python dir)
fasta = base / '..' / 'data' / '1_raw' / 'Homo_sapients.mtDNA.fasta'
fasta = fasta.resolve()
print('Using fasta:', fasta)

# run analyze
try:
    out = analyze_mtDNA_repeats(str(fasta), pos=8473, ref_nuc='T', alt_nuc='C', output_path=str(base / '..' / 'data' / '2_derived'))
    print('analyze returned:', out)
except Exception as e:
    print('analyze_mtDNA_repeats failed:', e)
    sys.exit(3)

# compare generated files with existing R-derived files if they exist
generated_all = out.get('all_repeats')
generated_summary = out.get('summary_top5_major_arc')

r_dir = (base / '..' / 'data' / '2_derived' / 'pos8473_TtoC').resolve()
r_all = r_dir / 'my_mtDNA_repeat_all_repeats.csv'
r_summary = r_dir / 'my_mtDNA_repeat_major_arc_summary_top5.csv'

print('Generated files:', generated_all, generated_summary)
print('Reference files:', r_all.exists(), r_summary.exists())

if r_all.exists() and Path(generated_all).exists():
    try:
        ga = pd.read_csv(generated_all)
        ra = pd.read_csv(r_all)
        print('all_repeats shapes gen/ref:', ga.shape, ra.shape)
        # quick compare: number of differing rows via merge anti-join
        merged = ga.merge(ra.drop_duplicates(), how='outer', indicator=True)
        diff_count = (merged['_merge'] != 'both').sum()
        print('Rows that differ (symmetric):', diff_count)
    except Exception as e:
        print('Error comparing all_repeats:', e)

if r_summary.exists() and Path(generated_summary).exists():
    try:
        gs = pd.read_csv(generated_summary)
        rs = pd.read_csv(r_summary)
        print('summary shapes gen/ref:', gs.shape, rs.shape)
        merged = gs.merge(rs.drop_duplicates(), how='outer', indicator=True)
        diff_count = (merged['_merge'] != 'both').sum()
        print('Summary rows that differ (symmetric):', diff_count)
    except Exception as e:
        print('Error comparing summary:', e)

print('Done')
