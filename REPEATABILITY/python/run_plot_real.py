#!/usr/bin/env python3
"""Load a repeats CSV and plot it using the viz helper. Saves a PNG.
"""
from pathlib import Path
import argparse
import pandas as pd
from src.viz import plot_circular_repeats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('csv', type=Path, help='Repeats CSV (mt_repeat_posXXX_all_repeats.csv)')
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parents[1] / 'viz_demo_out')
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    outdir = args.out
    outdir.mkdir(parents=True, exist_ok=True)
    plot_circular_repeats(df)
    print(f'Plotted {args.csv} to screen (and ensured {outdir} exists)')


if __name__ == '__main__':
    main()
