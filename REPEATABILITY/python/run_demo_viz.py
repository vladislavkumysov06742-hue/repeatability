import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.show = lambda *a, **k: None
import pandas as pd
import numpy as np
from pathlib import Path

from src import viz

outdir = Path('.').resolve() / 'viz_demo_out'
outdir.mkdir(exist_ok=True)

# Demo 1: bar chart
df = pd.DataFrame({'EffectiveLength':[10,12,15,18], 'Ref_Count':[5,3,2,1], 'Alt_Count':[1,2,3,4]})
viz.plot_repeatability_table(df)
plt.savefig(outdir / 'repeatability_bar.png')
plt.close()

# Demo 2: circular arcs
repeats = pd.DataFrame([
    {'pos':8473, 'repeat.start':100, 'repeat.end':200, 'EffectiveLength':18, 'RefAlt':'Ref'},
    {'pos':8473, 'repeat.start':300, 'repeat.end':500, 'EffectiveLength':12, 'RefAlt':'Alt'},
    {'pos':8473, 'repeat.start':1000, 'repeat.end':1100, 'EffectiveLength':16, 'RefAlt':'Alt'},
])
viz.plot_circular_repeats(repeats, radius=1.0, highlight_pos=8473)
plt.savefig(outdir / 'circular_repeats.png')
plt.close()

print('Saved demo plots to', outdir)
