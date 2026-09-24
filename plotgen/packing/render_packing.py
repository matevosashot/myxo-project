"""Render the packed capsules around a +1/2 defect from a saved shape file.

Usage:
    python render_packing.py [cells.json] [out.pdf]

The shape file holds, for each cell, the spine polyline (x, y nodes). A cell is the
set of points within width/2 of its spine. Only numpy and matplotlib are needed.
"""
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

# ---------------------------------------------------------------- appearance
FIGSIZE = 6.0                  # inches (square)
CELL_FACE = '#bcd8e8'
CELL_EDGE = '#1b3a5c'
CELL_LW = 1.0
FIXED_FACE = CELL_FACE         # the tail cell; change to highlight it
SHOW_STREAMLINES = True
STREAM_COLOR = '0.5'
STREAM_LW = 1.0
STREAM_P = np.concatenate([[0.0], np.arange(0.36, 6.0, 0.36)])   # streamline labels p
SHOW_BOX = True
BOX_COLOR = 'k'
BOX_LW = 2.0
SHOW_CORE = False
CORE_COLOR = 'red'
CORE_SIZE = 20
MARGIN = 0.02                  # axes margin around the box, in units of W
CAP_POINTS = 16                # points per half-circle cap
DPI = 300                      # for raster outputs (.png)


def capsule_polygon(spine, radius, nc=CAP_POINTS):
    """Outline of the set of points within `radius` of the spine polyline."""
    t = np.diff(spine, axis=0)
    t /= np.linalg.norm(t, axis=1, keepdims=True)
    tn = np.vstack([t[:1], 0.5 * (t[:-1] + t[1:]), t[-1:]])
    tn /= np.linalg.norm(tn, axis=1, keepdims=True)
    nrm = np.stack([-tn[:, 1], tn[:, 0]], 1)
    left, right = spine + radius * nrm, spine - radius * nrm
    a_end = np.arctan2(tn[-1, 1], tn[-1, 0])
    a_st = np.arctan2(tn[0, 1], tn[0, 0])
    ce = np.linspace(a_end + np.pi / 2, a_end - np.pi / 2, nc)
    cs = np.linspace(a_st - np.pi / 2, a_st - 3 * np.pi / 2, nc)
    cap_e = spine[-1] + radius * np.stack([np.cos(ce), np.sin(ce)], 1)
    cap_s = spine[0] + radius * np.stack([np.cos(cs), np.sin(cs)], 1)
    return np.vstack([left, cap_e, right[::-1], cap_s])


def streamline(p, extent):
    """Field line y^2 = 2 p x + p^2 (p = 0: the tail ray y = 0, x > 0)."""
    if p == 0:
        return np.array([[0.0, 0.0], [extent, 0.0]])
    y = np.linspace(-extent, extent, 2001)
    return np.stack([(y ** 2 - p ** 2) / (2 * p), y], 1)


def render(data, fname):
    W = data['W']
    radius = data['width'] / 2
    fig, ax = plt.subplots(figsize=(FIGSIZE, FIGSIZE))
    clip = Rectangle((-W, -W), 2 * W, 2 * W, transform=ax.transData)
    if SHOW_STREAMLINES:
        for p in STREAM_P:
            c = streamline(p, 3 * W)
            ln, = ax.plot(c[:, 0], c[:, 1], color=STREAM_COLOR, lw=STREAM_LW, zorder=1)
            ln.set_clip_path(clip)
    for cell in data['cells']:
        spine = np.asarray(cell['spine'], float)
        face = FIXED_FACE if cell['fixed'] else CELL_FACE
        pg = Polygon(capsule_polygon(spine, radius), closed=True,
                     fc=face, ec=CELL_EDGE, lw=CELL_LW, zorder=2)
        ax.add_patch(pg)
        pg.set_clip_path(clip)
    if SHOW_CORE:
        ax.scatter([0], [0], s=CORE_SIZE, color=CORE_COLOR, zorder=5)
    if SHOW_BOX:
        ax.add_patch(Rectangle((-W, -W), 2 * W, 2 * W, fc='none', ec=BOX_COLOR, lw=BOX_LW, zorder=4))
    lim = W * (1 + MARGIN)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.savefig(fname, bbox_inches='tight', dpi=DPI)
    plt.close(fig)


if __name__ == '__main__':
    src = sys.argv[1] if len(sys.argv) > 1 else 'cells_g0_polished_smooth.json'
    out = sys.argv[2] if len(sys.argv) > 2 else 'g0_polished_smooth.pdf'
    with open(src) as f:
        render(json.load(f), out)
    print('wrote', out)
