#!/usr/bin/env python3
"""Visualise the N dependence of the fluctuations_dirichlet.wls run time.

Reads the CSVs written by parse_log.py and produces

    fig1_stacked_absolute.png   wall time vs Nint, stacked by stage (+ small-N zoom)
    fig2_stacked_fraction.png   the same as a 100% stack -- where the time goes
    fig3_scaling_loglog.png     log-log scaling with fitted power laws

Only fbox=20 runs are plotted (the single fbox=10 run is a different physical
resolution and is not comparable); it stays in the CSVs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import (AutoMinorLocator, FixedLocator,
                               FuncFormatter, LogLocator)

HERE = Path(__file__).resolve().parent
FBOX = 20
ZOOM_MAX_N = 61          # second panel of fig1 zooms on Nint <= this

# The batch ran each resolution twice, at even and at odd Nint (odd puts the
# origin on the grid, even does not).  Only the even family is plotted; the odd
# runs stay in the CSVs.
DROP_ODD_NINT = True

# ---------------------------------------------------------------------------
# Stage grouping
# ---------------------------------------------------------------------------
# The two standalone residual evaluations are the same operation done twice;
# they are merged so the stack has one band per distinct kind of work.
MERGE = {
    "[6] residual eval #1 (unrefined)": "[6] residual evals (2x)",
    "[6] residual eval #3 (verification)": "[6] residual evals (2x)",
}

# Bands drawn explicitly = every stage whose share of tracked time reaches
# PEAK_SHARE_MIN in at least one run.  Everything else is folded into "other".
PEAK_SHARE_MIN = 2.5

# Fixed hue order -- a stage keeps its colour no matter its rank in a given run.
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                 "#e87ba4", "#008300", "#4a3aa7"]
OTHER_COLOR = "#8a8a85"

INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8a85"
GRID = "#e4e4e0"
GRID_MINOR = "#eae9e4"

SHORT = {
    "[6] Eigensystem (dense, nonsymmetric)": "[6] Eigensystem",
    "[6] refinement step (incl. residual eval #2)": "[6] refinement step",
    "[6] back-substitute": "[6] back-substitute",
    "[6] residual evals (2x)": "[6] residual evals",
    "[7] plot + export": "[7] plot + export",
    "[1] steady state (FEM, N-independent)": "[1] steady state (FEM)",
    "[3] kappa(S) inversion + g2,g3": "[3] kappa(S) inversion",
    "[7] save covariance (binary)": "[7] save covariance",
    "other": "other",
}

EIG = "[6] Eigensystem (dense, nonsymmetric)"
FEM = "[1] steady state (FEM, N-independent)"
PLOT = "[7] plot + export"
SAVE = "[7] save covariance (binary)"


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load():
    long = pd.read_csv(HERE / "timings_long.csv")
    long = long[long.fbox == FBOX].copy()
    if DROP_ODD_NINT:
        long = long[long.Nint % 2 == 0]
    long["stage"] = long["stage"].replace(MERGE)
    long = (long.groupby(["run_id", "Nint", "n", "stage"], as_index=False)
                .agg(seconds=("seconds", "sum")))
    total = long.groupby("Nint")["seconds"].transform("sum")
    long["share"] = 100 * long["seconds"] / total
    return long.sort_values(["Nint", "seconds"], ascending=[True, False])


def choose_bands(long):
    peak = long.groupby("stage")["share"].max().sort_values(ascending=False)
    bands = [s for s in peak.index if peak[s] >= PEAK_SHARE_MIN]
    return bands[:len(SERIES_COLORS)]


def pivot(long, bands):
    """Wide table Nint x band (seconds); everything else summed into 'other'."""
    w = long.pivot_table(index="Nint", columns="stage",
                         values="seconds", aggfunc="sum").fillna(0.0)
    other = w.drop(columns=bands).sum(axis=1)
    w = w[bands].copy()
    w["other"] = other
    return w


# ---------------------------------------------------------------------------
# styling helpers
# ---------------------------------------------------------------------------
def hours_fmt(v, _pos=None):
    """Decade tick label for an hours axis spanning 1e-5 h to 1e2 h."""
    if v <= 0:
        return ""
    exp = int(round(np.log10(v)))
    if abs(v - 10.0**exp) > 1e-9 * v:      # minor tick -- leave unlabelled
        return ""
    if -1 <= exp <= 2:
        return f"{10.0**exp:g}"
    return f"$10^{{{exp}}}$"


def style_axes(ax, *, grid_axis="y", minor_axis=None):
    ax.set_facecolor("white")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9, length=3, width=0.8)
    ax.tick_params(which="minor", colors=INK_MUTED, length=2, width=0.6)
    ax.grid(True, which="major", axis=grid_axis, color=GRID, lw=0.8, zorder=0)
    if minor_axis:
        ax.grid(True, which="minor", axis=minor_axis, color=GRID_MINOR,
                lw=0.6, zorder=0)
    ax.set_axisbelow(True)


def titles(fig, title, subtitle, x=0.01):
    """Figure-level title block, kept clear of the twin top axis."""
    return [
        fig.text(x, 0.99, title, fontsize=13, weight="bold", color=INK,
                 ha="left", va="top"),
        fig.text(x, 0.935, subtitle, fontsize=9.5, color=INK_2,
                 ha="left", va="top", linespacing=1.5),
    ]


def categorical_top_axis(ax, nints, *, label=True):
    """Twin top axis labelling the state dimension n = 3*Nint^2."""
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    top.set_xticks(range(len(nints)))
    rot = 0 if len(nints) <= 8 else 45
    top.set_xticklabels([f"{3 * N**2:,}" for N in nints],
                        rotation=rot, ha="center" if rot == 0 else "left",
                        fontsize=7.5, color=INK_MUTED)
    if label:
        top.set_xlabel("state dimension  n = 3 N$^2$", fontsize=8.5,
                       color=INK_MUTED, labelpad=6)
    for side in ("left", "right", "bottom"):
        top.spines[side].set_visible(False)
    top.spines["top"].set_color(GRID)
    top.tick_params(colors=INK_MUTED, length=3, width=0.8)
    return top


def legend_right(fig, ax, *, title="stage", reverse=True):
    handles, labels = ax.get_legend_handles_labels()
    if reverse:                       # stacks read bottom-up; match the bars
        handles, labels = handles[::-1], labels[::-1]
    leg = fig.legend(handles, labels, loc="upper left",
                     bbox_to_anchor=(1.005, 0.93), frameon=False, fontsize=9,
                     handlelength=1.1, handleheight=1.1, labelcolor=INK_2)
    leg.set_title(title, prop={"size": 9, "weight": "bold"})
    leg.get_title().set_color(INK)
    return leg


def save(fig, path, extra):
    fig.savefig(path, dpi=200, facecolor="white",
                bbox_inches="tight", bbox_extra_artists=extra)
    fig.savefig(path.with_suffix(".pdf"), facecolor="white",
                bbox_inches="tight", bbox_extra_artists=extra)
    plt.close(fig)
    print(f"wrote {path.name}")


# ---------------------------------------------------------------------------
# stacked bars
# ---------------------------------------------------------------------------
def draw_stack(ax, data, colors, *, fraction, label=True):
    """Draw one stacked-bar panel.  data: DataFrame indexed by Nint."""
    nints = list(data.index)
    x = np.arange(len(nints))
    plotted = 100 * data.div(data.sum(axis=1), axis=0) if fraction else data
    ceiling = 100.0 if fraction else plotted.sum(axis=1).max()

    bottom = np.zeros(len(nints))
    for col, color in zip(plotted.columns, colors):
        h = plotted[col].to_numpy()
        # 1.5pt surface gap between stacked segments
        ax.bar(x, h, bottom=bottom, width=0.72, color=color,
               edgecolor="white", linewidth=1.5,
               label=SHORT.get(col, col) if label else "_nolegend_", zorder=2)
        for xi, (hi, bi) in enumerate(zip(h, bottom)):
            if hi / ceiling > 0.085:      # selective direct labels only
                txt = f"{hi:.0f}%" if fraction else (
                    f"{hi:.0f}s" if hi >= 10 else f"{hi:.1f}s")
                ax.text(xi, bi + hi / 2, txt, ha="center", va="center",
                        fontsize=7.5, color="white", zorder=3)
        bottom += h

    ax.set_xticks(x)
    ax.set_xticklabels(nints, fontsize=9)
    ax.set_xlim(-0.7, len(nints) - 0.3)
    if fraction:
        ax.set_ylim(0, 100)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    style_axes(ax, minor_axis="y")
    return nints


def fig_absolute(w, bands, path):
    """Absolute seconds: full range plus a zoom on the small-N runs."""
    colors = SERIES_COLORS[:len(bands)] + [OTHER_COLOR]
    zoom = w[w.index <= ZOOM_MAX_N]

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(11.5, 5.4), gridspec_kw={"width_ratios": [1.25, 1]})
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(top=0.74, bottom=0.13, left=0.07, right=0.99, wspace=0.22)

    nA = draw_stack(axA, w[bands + ["other"]], colors, fraction=False)
    axA.set_ylabel("tracked wall time  [s]", fontsize=10, color=INK_2)
    axA.set_xlabel("interior grid points per axis, Nint", fontsize=10, color=INK_2)
    categorical_top_axis(axA, nA)
    axA.text(0.02, 0.97, "all runs", transform=axA.transAxes, fontsize=9.5,
             weight="bold", color=INK_2, va="top")

    nB = draw_stack(axB, zoom[bands + ["other"]], colors, fraction=False,
                    label=False)
    axB.set_xlabel("interior grid points per axis, Nint", fontsize=10, color=INK_2)
    categorical_top_axis(axB, nB)
    axB.text(0.02, 0.97, f"zoom: Nint $\\leq$ {zoom.index.max()}",
             transform=axB.transAxes, fontsize=9.5, weight="bold",
             color=INK_2, va="top")

    n_max, t_max = w.index.max(), w.sum(axis=1).max()
    txt = titles(fig,
                 "Run time vs grid resolution, by stage",
                 "fbox=20, Dirichlet BC, even Nint — tracked wall time per run.\n"
                 f"Left: the full range, which N={n_max} ({t_max:,.0f} s) "
                 "dominates. Right: the same stack for the runs it flattens.")
    leg = legend_right(fig, axA)
    save(fig, path, [leg, *txt])


def fig_fraction(w, bands, path):
    colors = SERIES_COLORS[:len(bands)] + [OTHER_COLOR]
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(top=0.74, bottom=0.13, left=0.09, right=0.99)

    nints = draw_stack(ax, w[bands + ["other"]], colors, fraction=True)
    ax.set_ylabel("share of tracked wall time  [%]", fontsize=10, color=INK_2)
    ax.set_xlabel("interior grid points per axis, Nint", fontsize=10, color=INK_2)
    categorical_top_axis(ax, nints)

    txt = titles(fig, "Where the time goes, by stage",
                 "The same runs as a 100% stack: N-independent overhead\n"
                 "dominates below N≈40, the dense eigensolve above it.")
    leg = legend_right(fig, ax)
    save(fig, path, [leg, *txt])


# ---------------------------------------------------------------------------
# log-log scaling
# ---------------------------------------------------------------------------
FIT_FROM_N = 60          # power laws are fitted on the asymptotic part only
PROJECT_TO_N = 200       # dotted extrapolation runs from the last run to here
HOUR = 3600.0


def fit_powerlaw(N, t, *, from_n=FIT_FROM_N):
    """Least-squares fit of t ~ C * N^p over the asymptotic part of the range.

    Returns (p, C).  Both are nan when there is nothing to fit.
    """
    N, t = np.asarray(N, float), np.asarray(t, float)
    keep = (t > 0) & (N >= from_n)
    if keep.sum() < 2:
        return np.nan, np.nan
    p, logC = np.polyfit(np.log(N[keep]), np.log(t[keep]), 1)
    return p, np.exp(logC)


def fit_slope(N, t, *, from_n=FIT_FROM_N):
    return fit_powerlaw(N, t, from_n=from_n)[0]


def fig_scaling(long, bands, w, path):
    """Log-log scaling in hours, with the fits extrapolated past the last run."""
    nints = np.array(w.index, float)
    n_last = nints.max()
    proj_N = np.geomspace(n_last, PROJECT_TO_N, 40)

    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(top=0.74, bottom=0.13, left=0.10, right=0.99)

    # save covariance is folded into "other" in the stacks but gets its own
    # curve here: it is the only I/O term and the clearest n^2 reference
    per_stage = long.pivot_table(index="Nint", columns="stage",
                                 values="seconds", aggfunc="sum").fillna(0.0)
    palette = dict(zip(bands, SERIES_COLORS))
    palette[SAVE] = OTHER_COLOR

    series = [("TOTAL tracked", w.sum(axis=1).to_numpy(), INK, 2.4, 6, 5)]
    series += [(SHORT.get(c, c), per_stage[c].to_numpy(),
                palette.get(c, OTHER_COLOR), 2.0, 5, 4)
               for c in list(bands) + [SAVE] if c in per_stage]

    for name, t_s, color, lw, ms, z in series:
        t_h = t_s / HOUR
        p, C = fit_powerlaw(nints, t_s)
        ax.plot(nints, t_h, color=color, lw=lw, marker="o", ms=ms, zorder=z,
                label=f"{name}  ($\\sim$N$^{{{p:.1f}}}$)")
        # dotted extrapolation of the same fit past the last completed run
        ax.plot(proj_N, C * proj_N**p / HOUR, color=color, lw=lw * 0.8,
                ls=":", zorder=z, label="_nolegend_")

    # the runs that did not finish: N=120 died in [6] densify A,Q at 27.8 GiB
    ax.axvspan(110, PROJECT_TO_N * 1.05, color="#f2f1ee", zorder=0)
    ax.axvline(110, color=INK_MUTED, lw=1, ls="--", zorder=1)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(18, PROJECT_TO_N * 1.05)
    ticks = [20, 30, 40, 60, 80, 100, 140, 200]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(
        FixedLocator([v for v in range(20, 201, 5) if v not in ticks]))
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_minor_locator(
        LogLocator(base=10.0, subs=tuple(range(2, 10)), numticks=200))
    ax.yaxis.set_major_formatter(FuncFormatter(hours_fmt))
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
    ax.set_xlabel("interior grid points per axis, Nint", fontsize=10, color=INK_2)
    ax.set_ylabel("wall time  [hours]", fontsize=10, color=INK_2)
    style_axes(ax, grid_axis="both", minor_axis="both")

    sec = ax.secondary_xaxis("top", functions=(lambda N: 3 * N**2,
                                               lambda n: np.sqrt(n / 3)))
    sec.xaxis.set_major_locator(FixedLocator([3 * N**2 for N in ticks]))
    sec.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(round(v)):,}"))
    sec.xaxis.set_minor_locator(FixedLocator([]))
    sec.set_xlabel("state dimension  n = 3 N$^2$", fontsize=8.5,
                   color=INK_MUTED, labelpad=6)
    sec.tick_params(colors=INK_MUTED, labelsize=8, length=3, width=0.8)
    sec.spines["top"].set_color(GRID)

    # annotate the projected total at the grid sizes that were attempted
    p_tot, C_tot = fit_powerlaw(nints, w.sum(axis=1).to_numpy())
    for N in (120, 140, 200):
        t_h = C_tot * N**p_tot / HOUR
        ax.plot([N], [t_h], marker="o", ms=5, mfc="white", mec=INK, mew=1.6,
                zorder=6)
        ax.annotate(f"N={N}: ~{t_h:,.0f} h", xy=(N, t_h),
                    xytext=(-6, 9), textcoords="offset points",
                    fontsize=8.5, color=INK, ha="right")
    ax.text(112, ax.get_ylim()[0] * 1.6,
            "N ≥ 120 OOM'd\n(27.8 GiB at densify)",
            fontsize=8.5, color=INK_2, ha="left", va="bottom")

    txt = titles(fig, "What the cost scales as — and what it would cost",
                 "Solid: measured. Dotted: each stage's own power law, fitted "
                 "over N ≥ 60 and extrapolated.\nA dense n×n eigensolve is "
                 "O(n³) = O(N⁶); writing the n×n covariance is O(n²) = O(N⁴).")
    leg = legend_right(fig, ax, title="stage (fitted slope)", reverse=False)
    save(fig, path, [leg, *txt])


# ---------------------------------------------------------------------------
def main():
    long = load()
    bands = choose_bands(long)
    w = pivot(long, bands)

    print(f"explicit bands (peak share >= {PEAK_SHARE_MIN}%): {len(bands)}")
    for b in bands:
        print(f"   {SHORT.get(b, b)}")

    fig_absolute(w, bands, HERE / "fig1_stacked_absolute.png")
    fig_fraction(w, bands, HERE / "fig2_stacked_fraction.png")
    fig_scaling(long, bands, w, HERE / "fig3_scaling_loglog.png")

    tot = w.sum(axis=1)
    print("\nNint    n        total[s]    eigen%   FEM+plot%")
    for N in w.index:
        print(f"{N:<6}  {3*N**2:<8} {tot[N]:>8.1f}  {100*w.loc[N, EIG]/tot[N]:>9.1f}  "
              f"{100*(w.loc[N, FEM]+w.loc[N, PLOT])/tot[N]:>10.1f}")
    print("\nfitted exponents (N >= 60):")
    per = long.pivot_table(index="Nint", columns="stage", values="seconds",
                           aggfunc="sum").fillna(0.0)
    Nv = np.array(w.index, float)
    for col in list(bands) + [SAVE]:
        print(f"   {SHORT.get(col, col):<24} N^{fit_slope(Nv, per[col].to_numpy()):.2f}")
    print(f"   {'TOTAL tracked':<24} N^{fit_slope(Nv, tot.to_numpy()):.2f}")


if __name__ == "__main__":
    main()
