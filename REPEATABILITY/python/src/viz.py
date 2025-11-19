import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from math import cos, sin, pi


def plot_repeatability_table(df: pd.DataFrame):
    if df.empty:
        print("No data to plot")
        return
    df.plot(x="EffectiveLength", y=["Ref_Count", "Alt_Count"], kind="bar")
    plt.tight_layout()
    plt.show()


def plot_circular_repeats(repeats_df: pd.DataFrame, radius: float = 1.0, highlight_pos: int = None):
    """Draw circular arcs similar to the R visualization script.

    repeats_df: expects columns pos, repeat.start, repeat.end, EffectiveLength, RefAlt
    """
    if repeats_df is None or repeats_df.empty:
        print("No repeats to plot")
        return

    def angle_for_pos(pos):
        MTDNA_LENGTH = 16568
        return (0.5 - (pos - 1) / MTDNA_LENGTH) * 2 * pi

    fig, ax = plt.subplots(figsize=(8, 8))
    # draw circle
    theta = np.linspace(0, 2 * pi, 400)
    ax.plot(np.cos(theta) * radius, np.sin(theta) * radius, color="gray", linewidth=2)

    # draw variant position marker (use first pos if not provided)
    pos_val = int(repeats_df['pos'].iloc[0]) if highlight_pos is None else int(highlight_pos)
    ang = angle_for_pos(pos_val)
    x_pos, y_pos = cos(ang) * radius, sin(ang) * radius
    ax.scatter([x_pos], [y_pos], color="goldenrod", s=80, zorder=5)
    ax.text(x_pos, y_pos + 0.08, str(pos_val), ha='center', va='bottom', fontweight='bold')

    # helper to create quadratic bezier points between two circle points
    def bezier_arc(x0, y0, x1, y1, curve=0.37, n=80):
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        dx, dy = x1 - x0, y1 - y0
        ctrl_x = xm - dy * curve
        ctrl_y = ym + dx * curve
        t = np.linspace(0, 1, n)
        x = (1 - t) ** 2 * x0 + 2 * (1 - t) * t * ctrl_x + t ** 2 * x1
        y = (1 - t) ** 2 * y0 + 2 * (1 - t) * t * ctrl_y + t ** 2 * y1
        return x, y

    # draw arcs
    for _, row in repeats_df.iterrows():
        rep_mid = (int(row['repeat.start']) + int(row['repeat.end'])) / 2.0
        ang_rep = angle_for_pos(rep_mid)
        x_rep, y_rep = cos(ang_rep) * radius, sin(ang_rep) * radius
        x_arc, y_arc = bezier_arc(x_pos, y_pos, x_rep, y_rep)
        color = '#3182bd' if row.get('RefAlt', 'Ref') == 'Ref' else '#de2d26'
        thickness = 7 if row.get('EffectiveLength', 0) >= 16 else 4
        ax.plot(x_arc, y_arc, color=color, linewidth=thickness, alpha=0.8)

    ax.set_aspect('equal')
    ax.axis('off')
    plt.show()
