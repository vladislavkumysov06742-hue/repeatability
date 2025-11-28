#!/usr/bin/env python3
"""Demo visualization using `src.viz.plot_circular_repeats`.
If a real repeats CSV is provided it will be used, otherwise a small
synthetic dataset is drawn.
"""
import argparse
from pathlib import Path
import pandas as pd
from src.viz import plot_circular_repeats


def make_demo_df():
    # small synthetic set of repeats for demo
    data = [
        {"pos": 8012, "repeat.start": 100, "repeat.end": 120, "EffectiveLength": 15, "RefAlt": "Ref"},
        {"pos": 8012, "repeat.start": 6000, "repeat.end": 6016, "EffectiveLength": 18, "RefAlt": "Alt"},
        {"pos": 8012, "repeat.start": 700, "repeat.end": 712, "EffectiveLength": 12, "RefAlt": "Ref"},
    ]
    return pd.DataFrame(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', type=Path, help='Path to repeats CSV to plot')
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parents[1] / 'viz_demo_out' / 'demo.png')
    args = parser.parse_args()

    if args.csv and args.csv.exists():
        df = pd.read_csv(args.csv)
    else:
        df = make_demo_df()

    plot_circular_repeats(df)
    print('Plotted demo (not saving by default)')


if __name__ == '__main__':
    main()
