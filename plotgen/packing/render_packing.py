"""Render the packed capsules around a +1/2 defect from a saved shape file.

Usage:
    python render_packing.py [cells.json] [out.pdf]
    python render_packing.py --out f.pdf --show-streams --pop-cell
    python render_packing.py --guides        # coordinate grid + p labels, to pick STREAM_AT

Run as a script, the feature flags are authoritative: anything not asked for is off,
so plain `--out f.pdf` gives bare cells whatever the SHOW_* constants below say. Those
constants set the appearance, and are the defaults when importing render() instead.

The shape file holds, for each cell, the spine polyline (x, y nodes). A cell is the
set of points within width/2 of its spine. Only numpy and matplotlib are needed.
"""
import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.colors import to_rgba
import matplotlib.patheffects as pe

# ---------------------------------------------------------------- appearance
FIGSIZE = 6.0                  # inches (square)
CELL_FACE = '#bcd8e8'
CELL_EDGE = '#1b3a5c'
CELL_LW = 0.5
CELL_ALPHA = 0.3              # fill only; outlines stay opaque
FIXED_FACE = CELL_FACE         # the tail cell; change to highlight it
SHOW_STREAMLINES = True
STREAM_COLOR = 'red'
STREAM_LW = 3
STREAM_ALPHA = 0.4
STREAM_ZORDER = 3              # 3: on top of the cells; 1: behind them
# One entry per streamline, in either form:
#   p        -- the line crosses the y axis at +-p and skirts the core at distance p/2
#   (x, y)   -- the line through that point
# p = 0 is the tail ray y = 0, x > 0. Only p < 2 W enters the box.
STREAM_AT = [0.0, 0.12, 0.36, 0.75, 1.2,  1.80]
SHOW_GUIDES = False            # True: axes ticks + p labels, to pick STREAM_AT by eye
SHOW_BOX = True
BOX_COLOR = 'k'
BOX_LW = 2.0
POP_CELLS = True               # lift chosen cells above the streamlines
POP_AT = [(-0.25, 0.4)]        # one point per cell to lift; the nearest cell wins
POP_FACE = '#9dd4f5'
POP_ALPHA = 1.0
POP_LW = 1.6
POP_ZORDER = 3.5               # above the streamlines (3), below the box (4)
SHOW_CORE = True               # red dot at the origin; drawn in every frame
CORE_COLOR = 'red'
CORE_SIZE = 70
SHOW_AXES_KEY = True           # x-y direction arrows in the bottom-left corner
KEY_POS = (-0.92, -0.92)       # where the two arrows meet, in units of W
KEY_LEN = 0.26                 # arrow length, in units of W
KEY_COLOR = 'k'
KEY_LW = 2.5
KEY_FS = 16
MARGIN = 0.005                  # axes margin around the box, in units of W
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


def stream_p(spec):
    """Label p of a STREAM_AT entry: a bare p, or a point (x, y) to pass through."""
    if np.ndim(spec) == 0:
        return float(spec)
    x, y = spec
    return float(np.hypot(x, y) - x)


def label_streamline(ax, c, p, W, i):
    """Tag a curve with its p, at its topmost point in the box; stagger to avoid overlap."""
    inside = c[(abs(c[:, 0]) <= W) & (abs(c[:, 1]) <= W)]
    if not len(inside):
        return
    ax.annotate(f'{p:.3g}', inside[np.argmax(inside[:, 1])], zorder=6,
                textcoords='offset points', xytext=(0, -3 - 10 * (i % 2)),
                color='crimson', fontsize=7, ha='center', va='top',
                bbox=dict(fc='w', ec='none', alpha=0.7, pad=0.5))


def pick_cells(cells, points):
    """Indices of the cells whose spines run nearest to each of `points`."""
    return {int(np.argmin([np.hypot(*(np.asarray(c['spine'], float) - q).T).min()
                           for c in cells])) for q in points}


def axes_key(ax, W):
    """Corner marker: arrows along +x and +y, with labels."""
    x0, y0 = KEY_POS[0] * W, KEY_POS[1] * W
    L = KEY_LEN * W
    halo = [pe.withStroke(linewidth=3, foreground='w')]
    for dx, dy, lab, off, ha, va in ((L, 0, 'x', (5, 0), 'left', 'center'),
                                     (0, L, 'y', (0, 5), 'center', 'bottom')):
        ax.annotate('', (x0 + dx, y0 + dy), xytext=(x0, y0), zorder=6,
                    arrowprops=dict(arrowstyle='-|>', color=KEY_COLOR, lw=KEY_LW,
                                    shrinkA=0, shrinkB=0, mutation_scale=14,
                                    # path_effects=halo
                                    ))
        ax.annotate(lab, (x0 + dx, y0 + dy), textcoords='offset points', xytext=off,
                    ha=ha, va=va, color=KEY_COLOR, fontsize=KEY_FS, style='italic',
                    zorder=6, path_effects=halo)


def render(data, fname):
    W = data['W']
    radius = data['width'] / 2
    fig, ax = plt.subplots(figsize=(FIGSIZE, FIGSIZE))
    clip = Rectangle((-W, -W), 2 * W, 2 * W, transform=ax.transData)
    if SHOW_STREAMLINES:
        for i, spec in enumerate(STREAM_AT):
            p = stream_p(spec)
            c = streamline(p, 3 * W)
            ln, = ax.plot(c[:, 0], c[:, 1], color=STREAM_COLOR, lw=STREAM_LW,
                          alpha=STREAM_ALPHA, solid_capstyle='round',
                          zorder=STREAM_ZORDER)
            ln.set_clip_path(clip)
            if SHOW_GUIDES:
                label_streamline(ax, c, p, W, i)
    popped = pick_cells(data['cells'], POP_AT) if POP_CELLS else set()
    for i, cell in enumerate(data['cells']):
        spine = np.asarray(cell['spine'], float)
        if i in popped:
            face, alpha, lw, z = POP_FACE, POP_ALPHA, POP_LW, POP_ZORDER
        else:
            face = FIXED_FACE if cell['fixed'] else CELL_FACE
            alpha, lw, z = CELL_ALPHA, CELL_LW, 2
        pg = Polygon(capsule_polygon(spine, radius), closed=True,
                     fc=to_rgba(face, alpha), ec=CELL_EDGE, lw=lw, zorder=z)
        ax.add_patch(pg)
        pg.set_clip_path(clip)
    if SHOW_CORE:
        ax.scatter([0], [0], s=CORE_SIZE, color=CORE_COLOR, zorder=5)
    if SHOW_BOX:
        ax.add_patch(Rectangle((-W, -W), 2 * W, 2 * W, fc='none', ec=BOX_COLOR, lw=BOX_LW, zorder=4))
    if SHOW_AXES_KEY:
        axes_key(ax, W)
    lim = W * (1 + MARGIN)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    if SHOW_GUIDES:
        for setter in (ax.set_xticks, ax.set_yticks):
            setter(np.linspace(-W, W, 9))
            setter(np.linspace(-W, W, 41), minor=True)
        ax.tick_params(labelsize=7, color='crimson', labelcolor='crimson')
        ax.tick_params(which='minor', length=2, color='crimson')
        ax.grid(color='crimson', alpha=0.3, lw=0.6)
        ax.grid(which='minor', color='crimson', alpha=0.15, lw=0.4)
        ax.set_axisbelow(False)
    else:
        ax.axis('off')
    fig.savefig(fname, bbox_inches='tight', dpi=DPI)
    plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('cells', nargs='?', default='cells_g0_polished_smooth.json')
    ap.add_argument('out_pos', nargs='?', metavar='out', help='output file')
    ap.add_argument('-o', '--out', dest='out', help='output file (same as the 2nd argument)')
    ap.add_argument('--show-streams', action='store_true', help='draw the field lines')
    ap.add_argument('--pop-cell', action='store_true', help='lift the POP_AT cells on top')
    ap.add_argument('--guides', action='store_true', help='coordinate grid + p labels')
    a = ap.parse_args()

    SHOW_STREAMLINES, POP_CELLS, SHOW_GUIDES = a.show_streams, a.pop_cell, a.guides
    if a.guides:
        SHOW_STREAMLINES, FIGSIZE = True, 10.0
    out = a.out or a.out_pos or ('guides.png' if a.guides else 'g0_polished_smooth.pdf')
    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    with open(a.cells) as f:
        render(json.load(f), out)
    print('wrote', out)
