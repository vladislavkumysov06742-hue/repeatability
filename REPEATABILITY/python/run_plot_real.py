import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.show = lambda *a, **k: None
import pandas as pd
from pathlib import Path
from src import viz

csv_path = Path('.').resolve() / '.pytest_tmp' / 'test_r_and_python_equivalence0' / 'pos8473_TtoC' / 'my_mtDNA_repeat_all_repeats.csv'
print('Looking for CSV at', csv_path)
if not csv_path.exists():
    print('CSV not found:', csv_path)
    raise SystemExit(1)

df = pd.read_csv(csv_path, dtype=str)
# Ensure numeric columns
for col in ['repeat.start','repeat.end','EffectiveLength','pos']:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

print('Loaded repeats:', len(df))
if df.empty:
    print('No repeats to plot')
    raise SystemExit(0)

viz.plot_circular_repeats(df, radius=1.0, highlight_pos=int(df['pos'].iloc[0]))
plt.savefig(Path('.').resolve() / 'viz_demo_out' / 'pos8473_real_repeats.png')
plt.close()
print('Saved viz to', Path('.').resolve() / 'viz_demo_out' / 'pos8473_real_repeats.png')
