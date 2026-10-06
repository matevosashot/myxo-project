#!/usr/bin/env python3
"""Nucleation barrier under a fluctuating pressure.

    python3 remote_kernel/plot_fluctuating_barier.py
    python3 remote_kernel/plot_fluctuating_barier.py --out /tmp/b.png --n 200

Writes a PDF by default.  The two bundles -- 2 x n translucent curves, the one
thing in the figure that makes a vector file big and slow to open -- are
rasterized at --raster-dpi; frame, ticks, text and the critical curve stay vector.

The style follows ./figures: a 176 pt wide canvas, 6 pt text, 1.44 pt curves,
an inward-tick frame, no title, and the palette those figures already use.

Port of the Mathematica snippet:

    params = {Ws -> 2000, Wc -> 0, Gamma -> 65000, h -> 0.5, Ac -> 3};
    e[A_, p_] = (Ws - Wc - p) A + (Gamma + Wc) h Sqrt[4 Pi A];

The free energy of a circular nucleus of area A, at pressure p: a bulk term
linear in A and a line-tension term going as the perimeter, 2 Sqrt(Pi A), times
the layer height h.  For p above the critical pressure the bulk term wins at
large A and the curve has a maximum -- the nucleation barrier.

WHERE pc COMES FROM -- the snippet leaves it undefined
------------------------------------------------------
pc is fixed by Ac: it is the pressure whose barrier top sits at A = Ac.
With b = (Gamma + Wc) h 2 Sqrt(Pi), de/dA = (Ws - Wc - p) + b / (2 Sqrt(A)),
so the top is at Sqrt(A*) = b / (2 (p - Ws + Wc)) and demanding A* = Ac gives

    pc = Ws - Wc + (Gamma + Wc) h Sqrt(Pi / Ac)        (35.258 kPa um here)

Override it with --pc if the intended definition is a different one.

WHAT IT SHOWS
-------------
Two bundles of curves around the critical one, from a pressure that fluctuates
by a few percent: p = pc / (1 + |x|) with x ~ N(0, 0.1) below, and
p = pc (1 + |x|) with x ~ N(0, 0.2) above.  The asymmetry is the point -- the
barrier top moves far more, and far faster, on the low-pressure side, so a
symmetric jitter in p produces a strongly skewed spread in barrier height.
"""

import argparse
import math
import os

# House style, measured off figures/final.pdf, ln_curve.pdf, strength.pdf:
# a 176 pt wide Mathematica export, every glyph at 6 pt (4 pt subscripts),
# data curves at 1.44 pt, a full inward-tick frame, no title, no grid, and the
# ColorData[97] palette.  Sizes below are in points on that 176 pt canvas, so
# the figure drops into the manuscript at the same scale as its siblings.
BLUE, GOLD, GREEN, RED = "#5e81b5", "#e19c24", "#8fb131", "#ff0000"
INK, FAINT, SURF = "#000000", "#cccccc", "#ffffff"

WIDTH_PT, HEIGHT_PT = 126.0, 126.0   # figures/final.pdf page size
LEGEND_XY = (0.99, 0.99)             # top-right corner, in axes fraction
LEGEND_PAD = (0.3, 0.3)              # x, y padding inside the legend plate
LEGEND_ALPHA = 0.1              # opacity of that white plate
LEGEND_EDGE_LW = 0.5                 # its border, thinner than the frame
FS = 6.0                             # every label in the house figures
LW_CURVE, LW_FRAME, LW_FAINT = 1.044, 0.8, 0.144

# manuscript.tex, tab:pcrit: energies in kPa um (= mN/m), lengths in um, so
# e(A, p) comes out in kPa um^3 for A in um^2
PARAMS = dict(Ws=2.0, Wc=0.0, gamma=65.0, h=0.5, Ac=3.0)


def energy(A, p, P):
    """e(A, p) = (Ws - Wc - p) A + (gamma + Wc) h sqrt(4 pi A)."""
    import numpy as np
    return ((P["Ws"] - P["Wc"] - p) * A
            + (P["gamma"] + P["Wc"]) * P["h"] * np.sqrt(4.0 * np.pi * A))


def critical_pressure(P):
    """The p whose barrier top sits at A = Ac."""
    return P["Ws"] - P["Wc"] + (P["gamma"] + P["Wc"]) * P["h"] * math.sqrt(math.pi / P["Ac"])


def barrier_top(p, P):
    """(A*, e(A*)) of the barrier at pressure p; None if the curve has no maximum."""
    c = P["Ws"] - P["Wc"] - p
    if c >= 0:
        return None
    b = (P["gamma"] + P["Wc"]) * P["h"] * 2.0 * math.sqrt(math.pi)
    A = (b / (-2.0 * c)) ** 2
    return A, energy(A, p, P)


def plot(out_path, P, pc, n, amax, ylim, sigma_lo, sigma_hi, seed, raster_dpi,
         legend_xy=LEGEND_XY, legend_pad=LEGEND_PAD, legend_alpha=LEGEND_ALPHA,
         width_pt=WIDTH_PT, height_pt=HEIGHT_PT):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import logging
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.ticker import MultipleLocator

    # CMU Serif registers at weight 500; the lookup warning is noise
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

    # Computer Modern, the face the Mathematica exports use
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["CMU Serif", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
        "font.size": FS,
    })

    rng = np.random.default_rng(seed)
    p_lo = pc / (1.0 + np.abs(rng.normal(0.0, sigma_lo, n)))   # below pc
    p_hi = pc * (1.0 + np.abs(rng.normal(0.0, sigma_hi, n//2)))   # above pc

    # p_lo = pc / (1.0 + np.abs(rng.exponential(sigma_lo/2, n)))   # below pc
    # p_hi = pc * (1.0 + np.abs(rng.exponential(sigma_hi/2, n)))   # above pc

    A = np.linspace(0.0, amax, 1200)
    def e(p):                                   # kPa um^3
        return energy(A, p, P)

    fig = plt.figure(figsize=(width_pt / 72.0, height_pt / 72.0))
    fig.patch.set_facecolor(SURF)
    rect = [0.195, 0.180, 0.745, 0.745]
    ax = fig.add_axes(rect)
    ax.set_facecolor(SURF)

    # The two bundles: one thin translucent curve per sampled pressure.  Each
    # bundle goes in as ONE rasterized LineCollection, so the 2n overlapping
    # translucent strokes -- the only thing here that would bloat a vector file
    # -- become two flat images, while the frame, ticks, text and p_c curve
    # around them stay vector.  (Note this is not ax.set_rasterization_zorder:
    # that would drag anything drawn below the data into the raster too.)
    for ps, colour in ((p_lo, BLUE), (p_hi, GOLD)):
        segs = [np.column_stack([A, e(p)]) for p in ps]
        ax.add_collection(LineCollection(segs, colors=colour, linewidths=0.5,
                                         alpha=0.05, rasterized=True, zorder=2))

    ax.axhline(0.0, color=FAINT, lw=LW_FAINT, zorder=1)
    ax.plot(A, e(pc), color=INK, lw=LW_CURVE, solid_capstyle="round", zorder=4)

    # the critical barrier top, at A = Ac by construction -- marked with the
    # faint reference line figures/ln_curve.pdf uses for the same job
    top = barrier_top(pc, P)
    if top:
        ax.axvline(top[0], color=FAINT, lw=LW_FAINT, zorder=1)
        # in the clear band below the bundles, beside the reference line
        tr = mpl.transforms.blended_transform_factory(ax.transData, ax.transAxes)
        ax.text(top[0], 0.40, r"  $A_{\rm cell}$", transform=tr, color=INK,
                fontsize=FS, ha="left", va="center", zorder=6)

    handles = [plt.Line2D([], [], color=BLUE, lw=LW_CURVE),
               plt.Line2D([], [], color=INK, lw=LW_CURVE),
               plt.Line2D([], [], color=GOLD, lw=LW_CURVE)]
    # Anchored by its top-right corner at legend_xy in axes fraction, on a
    # rounded translucent white plate so the bundles read through it.  The
    # boxstyle pad is drawn OUTSIDE the legend's own bbox, so back it out of
    # the anchor -- otherwise the plate spills past the frame.  legend_xy of
    # (1, 1) then means flush with the frame's top-right corner.
    bpad = legend_pad[0] * FS                      # points
    anchor = (legend_xy[0] - bpad / (width_pt * rect[2]),
              legend_xy[1] - bpad / (height_pt * rect[3]))
    leg = ax.legend(handles,
                    [r"$p<p_{\rm crit}$", r"$p=p_{\rm crit}$", r"$p>p_{\rm crit}$"],
                    loc="upper right", bbox_to_anchor=anchor,
                    bbox_transform=ax.transAxes,
                    frameon=True, fancybox=True,
                    facecolor=SURF, edgecolor=INK, fontsize=FS,
                    handlelength=1.5, handletextpad=0.5, labelspacing=0.35,
                    borderpad=legend_pad[1], borderaxespad=0.0)
    leg.set_zorder(6)
    # framealpha (and its rcParam fallback) sets ONE alpha on the whole patch,
    # which both overrides an RGBA facecolor and fades the border with it.  So
    # clear the patch alpha and carry the opacity in the face colour alone --
    # the border then stays solid at any --legend-alpha.
    fr = leg.get_frame()
    fr.set_boxstyle("round", pad=legend_pad[0], rounding_size=0.35)
    fr.set_alpha(None)
    fr.set_facecolor(mpl.colors.to_rgba(SURF, legend_alpha))
    fr.set_edgecolor(INK)
    fr.set_linewidth(LEGEND_EDGE_LW)

    ax.set_xlim(-1, amax)
    ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(MultipleLocator(3.0))
    ax.set_xlabel(r"Patch area  $A$  ($\mu\mathrm{m}^2$)", fontsize=FS, labelpad=3.5)
    ax.set_ylabel(r"$\Delta E(A)$  ($\mathrm{kPa}\cdot\mu\mathrm{m}^3$)", fontsize=FS, labelpad=3.5)

    # full frame; ticks inward on the left and bottom only
    for sp in ax.spines.values():
        sp.set_color(INK)
        sp.set_linewidth(LW_FRAME)
    ax.tick_params(which="major", direction="in", top=False, right=False,
                   length=2.6, width=0.6, color=INK, labelsize=FS, pad=3.0)

    fig.savefig(out_path, facecolor=SURF, dpi=raster_dpi)
    return top


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--out", default=os.path.join(here, "fluctuating_barrier.pdf"))
    ap.add_argument("--raster-dpi", type=float, default=800.0,
                    help="resolution of the rasterized bundle layer")
    ap.add_argument("--n", type=int, default=50, help="curves per bundle")
    ap.add_argument("--amax", type=float, default=15.0, help="right edge of the A axis")
    ap.add_argument("--ylim", type=float, nargs=2, default=(-25.0, 115.0),
                    metavar=("LO", "HI"), help="y range, in kPa um^3")
    ap.add_argument("--sigma-lo", type=float, default=0.1)
    ap.add_argument("--sigma-hi", type=float, default=0.05)
    ap.add_argument("--pc", type=float, default=None, help="override the critical pressure")
    ap.add_argument("--seed", type=int, default=1, help="for a reproducible sample")
    ap.add_argument("--legend-xy", type=float, nargs=2, default=LEGEND_XY,
                    metavar=("X", "Y"),
                    help="legend top-right corner, in axes fraction")
    ap.add_argument("--legend-pad", type=float, nargs=2, default=LEGEND_PAD,
                    metavar=("PX", "PY"),
                    help="x and y padding inside the legend plate, in font units")
    ap.add_argument("--legend-alpha", type=float, default=LEGEND_ALPHA,
                    help="opacity of the legend's white plate (0 = invisible)")
    ap.add_argument("--width-pt", type=float, default=WIDTH_PT)
    ap.add_argument("--height-pt", type=float, default=HEIGHT_PT)
    for k, v in PARAMS.items():
        ap.add_argument(f"--{k}", type=float, default=v)
    a = ap.parse_args()

    P = {k: getattr(a, k) for k in PARAMS}
    pc = a.pc if a.pc is not None else critical_pressure(P)
    top = plot(a.out, P, pc, a.n, a.amax, tuple(a.ylim), a.sigma_lo, a.sigma_hi, a.seed, a.raster_dpi,
               tuple(a.legend_xy), tuple(a.legend_pad), a.legend_alpha,
               a.width_pt, a.height_pt)
    print(f"pc = {pc:.4f}")
    if top:
        print(f"barrier top at A* = {top[0]:.4f}, e = {top[1]:.4f}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
