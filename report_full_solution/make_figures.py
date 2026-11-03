#!/usr/bin/env python3
"""Figures of the full-model report (report_full_solution/report.tex).

Data:   data/<tag>_m<k>.json and data/<tag>_runs.json, written by export_fields.wls
        from ../comoving_defect_full/out/<tag>.mx.
        Milestones of path A: m1 paper (all weights 0), m2 {wT}, m3 {wT, wR},
        m4 full {wT, wR, wH, wE}.  Path B: m1 paper, m2 {wR}, m3 {wR, wH, wE}, m4 full.
Output: figures/*.pdf.  House style of
        claude_experiments/drift_closure_2026_10_09/make_figure.py.
Run:    ../.venv/bin/python -I make_figures.py [--A A1] [--B B1]
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FIG = os.path.join(HERE, "figures")

plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{upgreek}\usepackage{amsmath}",
    "font.family": "serif",
    "font.size": 9,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": False, "ytick.right": False,   # ticks on the left and bottom axes only
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
})
BLUE = "#3f6fb0"      # full model (house blue)
ORANGE = "#d9600b"    # intermediate models
RED = "#b8322a"
GREEN = "#2e8540"
GRAY = "0.45"
UM = r"$\upmu$m"


def load(tag, k):
    path = os.path.join(DATA, f"{tag}_m{k}.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def load_runs(tag):
    path = os.path.join(DATA, f"{tag}_runs.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def pade_S(r, S0=0.9, ld=1.0):
    x = np.asarray(r)/ld
    return S0*x*np.sqrt((0.34 + 0.07*x**2)/(1 + 0.41*x**2 + 0.07*x**4))


def panel_label(ax, s, x=-0.17):
    ax.text(x, 1.02, r"\textbf{(%s)}" % s, transform=ax.transAxes, ha="left", va="bottom")


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    out = os.path.join(FIG, name)
    fig.savefig(out)
    plt.close(fig)
    print("wrote", os.path.relpath(out, HERE))


# ---------------------------------------------------------------- u along the paths
def fig_u_path(tagA, tagB):
    rA, rB = load_runs(tagA), load_runs(tagB)
    if rA is None:
        return
    zt = rA["zeta"]
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.6), sharey=True,
                            gridspec_kw=dict(wspace=0.12, left=0.085, right=0.985,
                                             bottom=0.13, top=0.9))

    def coordA(p):            # p = {zeta, wT, wR, wH, wE}
        z, wT, wR, wH, wE = p
        if abs(z - zt) > 1e-9:
            return None
        if wR == 0 and wH == 0:
            return wT
        if wH == 0:
            return 1 + wR
        return 2 + wH

    def coordB(p):            # paper -> wR -> (wH, wE) -> wT
        z, wT, wR, wH, wE = p
        if abs(z - zt) > 1e-9:
            return None
        if wH == 0 and wT == 0:
            return wR
        if wT == 0:
            return 1 + wH
        return 2 + wT

    for ax, run, coord, ticks, names in (
            (axs[0], rA, coordA, [0, 1, 2, 3],
             [r"paper", r"$+w_T$", r"$+w_R$", r"full"]),
            (axs[1], rB, coordB, [0, 1, 2, 3],
             [r"paper", r"$+w_R$", r"$+w_H,w_E$", r"full"])):
        if run is None:
            ax.set_visible(False)
            continue
        pts = [(coord(t["p"]), t["u"][0]) for t in run["trace"]]
        pts = sorted([p for p in pts if p[0] is not None])
        s, u = np.array(pts).T
        ax.axhline(0, color=GRAY, lw=0.6, zorder=0)
        ax.plot(s, u, "-", color=BLUE, lw=1.0, zorder=2)
        ax.plot(s, u, "o", ms=2.6, color=BLUE, zorder=3)
        for m in run["milestones"]:
            c = coord([zt] + m["w"])
            ax.plot([c], [m["u"][0]], "o", ms=5.5, mfc="white", mec="black", mew=0.9,
                    zorder=4)
        for t in ticks[1:-1]:
            ax.axvline(t, color=GRAY, lw=0.5, ls=":", zorder=0)
        ax.set_xticks(ticks)
        ax.set_xticklabels(names)
        ax.set_xlim(ticks[0] - 0.12, ticks[-1] + 0.12)
    axs[0].set_ylabel(r"$u_x$ (%s/min)" % UM)
    panel_label(axs[0], "a", x=0.0)
    panel_label(axs[1], "b", x=0.0)
    axs[0].text(0.03, 0.05, "path A", transform=axs[0].transAxes, color=GRAY)
    axs[1].text(0.03, 0.05, "path B", transform=axs[1].transAxes, color=GRAY)
    save(fig, "u_path.pdf")


# ---------------------------------------------------------------- S along the axes
def fig_S_axes(tagA):
    paper, full = load(tagA, 1), load(tagA, 4)
    mid = load(tagA, 3)
    if paper is None or full is None:
        return
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.6),
                            gridspec_kw=dict(wspace=0.28, left=0.085, right=0.985,
                                             bottom=0.17, top=0.9))
    for ax, rmax in ((axs[0], 5.0), (axs[1], None)):
        for d, ls in ((paper, "--"), (full, "-")):
            x = np.array(d["axisX"]["x"]); S = np.array(d["axisX"]["S"])
            y = np.array(d["axisY"]["y"]); Sy = np.array(d["axisY"]["S"])
            col = BLUE if d is full else "black"
            ax.plot(x[x >= 0], S[x >= 0], ls, color=col, lw=1.0)            # tail, +x
            ax.plot(-x[x <= 0], S[x <= 0], ls, color=col, lw=1.0, alpha=0.55)  # head, -x
            ax.plot(y[y >= 0], Sy[y >= 0], ls, color=col, lw=0.8, alpha=0.3)  # side, +y
        r = np.linspace(0, 20, 400)
        ax.plot(r, pade_S(r), ":", color=GRAY, lw=1.0)
        ax.set_xlabel(r"$r$ (%s)" % UM, labelpad=1)
        ax.set_ylabel(r"$S$")
        if rmax:
            ax.set_xlim(0, rmax)
            ax.set_ylim(0, 0.95)
        else:
            ax.set_xlim(0, 20)
            ax.set_ylim(0.78, 0.92)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color="black", ls="--", lw=1.0, label="paper"),
         Line2D([], [], color=BLUE, ls="-", lw=1.0, label="full"),
         Line2D([], [], color=GRAY, ls=":", lw=1.0, label=r"Pad\'e"),
         Line2D([], [], color="k", lw=1.0, alpha=1.0, label=r"tail ($+x$)"),
         Line2D([], [], color="k", lw=1.0, alpha=0.55, label=r"head ($-x$)"),
         Line2D([], [], color="k", lw=0.8, alpha=0.3, label=r"side ($+y$)")]
    axs[0].legend(handles=h, loc="lower right", frameon=False, ncol=2, handlelength=1.6,
                  columnspacing=0.8)
    panel_label(axs[0], "a", x=0.0)
    panel_label(axs[1], "b", x=0.0)
    save(fig, "S_axes.pdf")


# ---------------------------------------------------------------- director deviation
def fig_director_dev(tagA):
    paper, full = load(tagA, 1), load(tagA, 4)
    if paper is None or full is None:
        return
    fig, ax = plt.subplots(figsize=(3.375, 2.5))
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.18, top=0.96)
    cols = {0.5: GREEN, 1.0: BLUE, 2.0: ORANGE, 5.0: RED}
    for d, ls in ((paper, "--"), (full, "-")):
        for c in d["circles"]:
            r = c["r"]
            if r not in cols:
                continue
            ax.plot(np.degrees(c["phi"]), np.degrees(c["dev"]), ls, color=cols[r], lw=1.0,
                    label=(r"$r=%g$ %s" % (r, UM)) if d is full else None)
    ax.axhline(0, color=GRAY, lw=0.5)
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xlabel(r"$\phi$ (deg)", labelpad=1)
    ax.set_ylabel(r"$\theta-\phi/2$ (deg)")
    ax.legend(loc="upper left", frameon=False, handlelength=1.6)
    save(fig, "director_dev.pdf")


# ---------------------------------------------------------------- density on the axis
def fig_rho_axis(tagA):
    paper, full, mid = load(tagA, 1), load(tagA, 4), load(tagA, 2)   # mid: +w_T
    if paper is None or full is None:
        return
    fig, ax = plt.subplots(figsize=(3.375, 2.5))
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.18, top=0.96)
    R = paper["box"]
    x = np.linspace(-R, R, 801)
    ax.plot(x, -(pade_S(np.abs(x)) - np.abs(x)/R*pade_S(R))*np.sign(x), ":", color=GRAY,
            lw=1.0, label=r"closed form, Pad\'e")
    ax.plot(paper["axisX"]["x"], paper["axisX"]["rhoScaled"], "--", color="black", lw=1.0,
            label="paper")
    if mid is not None:
        ax.plot(mid["axisX"]["x"], mid["axisX"]["rhoScaled"], "-", color=ORANGE, lw=0.9,
                label=r"$+w_T$")
    ax.plot(full["axisX"]["x"], full["axisX"]["rhoScaled"], "-", color=BLUE, lw=1.1,
            label="full")
    ax.axhline(0, color=GRAY, lw=0.5)
    ax.set_xlim(-R, R)
    ax.set_xlabel(r"$x$ (%s), $y=0$" % UM, labelpad=1)
    ax.set_ylabel(r"$(\rho-1)\,B/\zeta$")
    ax.legend(loc="upper right", frameon=False, handlelength=1.6)
    save(fig, "rho_axis.pdf")


# ---------------------------------------------------------------- forces on the axis
def fig_forces_axis(tagA):
    paper, full = load(tagA, 1), load(tagA, 4)
    if paper is None or full is None:
        return
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.6), sharey=True,
                            gridspec_kw=dict(wspace=0.08, left=0.085, right=0.985,
                                             bottom=0.17, top=0.9))
    for ax, d, name in ((axs[0], paper, "paper"), (axs[1], full, "full")):
        x = np.array(d["axisX"]["x"])
        fa = np.array(d["axisX"]["fActive"]); fe = np.array(d["axisX"]["fElastic"])
        fp = np.array(d["axisX"]["fPressure"])
        el_on = name == "full"
        ax.plot(x, fa, "-", color=RED, lw=1.0, label=r"active $f^a_x$")
        if el_on:
            ax.plot(x, fe, "-", color=GREEN, lw=1.0, label=r"elastic $(\mathbf Q\!:\!\partial_x\mathbf H)$")
        ax.plot(x, fp, "-", color=BLUE, lw=1.0, label=r"pressure $-B\partial_x\rho$")
        tot = fa + fp + (fe if el_on else 0)
        ax.plot(x, tot, "--", color="black", lw=0.9, label=r"sum $=\xi_0 v_x$")
        ax.axhline(0, color=GRAY, lw=0.5)
        ax.set_xlim(-5, 5)
        ax.set_xlabel(r"$x$ (%s), $y=0$" % UM, labelpad=1)
        ax.text(0.97, 0.95, name, transform=ax.transAxes, ha="right", va="top", color=GRAY)
    axs[0].set_ylabel(r"force density (Pa)")
    axs[1].legend(loc="lower left", frameon=False, handlelength=1.6, fontsize=7.5)
    panel_label(axs[0], "a", x=0.0)
    panel_label(axs[1], "b", x=0.0)
    save(fig, "forces_axis.pdf")


# ---------------------------------------------------------------- maps near the core
def fig_maps(tagA):
    paper, full = load(tagA, 1), load(tagA, 4)
    if paper is None or full is None:
        return
    fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.35),
                            gridspec_kw=dict(wspace=0.45, left=0.06, right=0.95,
                                             bottom=0.17, top=0.93))
    mf, mp = full["maps"], paper["maps"]
    x, y = np.array(mf["x"]), np.array(mf["y"])
    ext = [x[0], x[-1], y[0], y[-1]]
    S = np.array(mf["S"]); dS = S - np.array(mp["S"]); H = np.array(mf["Habs"])
    for ax, Z, cmap, lab, sym in ((axs[0], S, "viridis", r"$S$ (full)", False),
                                  (axs[1], dS, "RdBu_r", r"$S_\mathrm{full}-S_\mathrm{paper}$", True),
                                  (axs[2], H, "magma", r"$|\mathbf H|$ (full)", False)):
        vmax = np.max(np.abs(Z)) if sym else None
        im = ax.imshow(Z, origin="lower", extent=ext, cmap=cmap,
                       vmin=(-vmax if sym else None), vmax=vmax)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=7)
        ax.set_title(lab, fontsize=9, pad=3)
        ax.set_xlabel(r"$x$ (%s)" % UM, labelpad=1)
        ax.set_xticks([-4, 0, 4]); ax.set_yticks([-4, 0, 4])
    axs[0].set_ylabel(r"$y$ (%s)" % UM)
    save(fig, "maps.pdf")


# ---------------------------------------------------------------- director lines
def director_lines(d, seeds, box=9.5, h=0.05, nmax=600, rcore=0.25):
    xs, ys = np.array(d["director"]["x"]), np.array(d["director"]["y"])
    Q1 = RegularGridInterpolator((ys, xs), np.array(d["director"]["Q1"]))
    Q2 = RegularGridInterpolator((ys, xs), np.array(d["director"]["Q2"]))

    def n_at(p):
        q1, q2 = Q1((p[1], p[0])), Q2((p[1], p[0]))
        th = 0.5*np.arctan2(q2, q1)
        return np.array([np.cos(th), np.sin(th)])

    lines = []
    for s in seeds:
        for sgn in (1, -1):
            p = np.array(s, float); prev = None; pts = [p.copy()]
            for _ in range(nmax):
                n1 = n_at(p)
                if prev is None:
                    n1 = sgn*n1
                elif n1 @ prev < 0:
                    n1 = -n1
                pm = p + 0.5*h*n1
                if np.max(np.abs(pm)) > box:
                    break
                n2 = n_at(pm)
                if n2 @ n1 < 0:
                    n2 = -n2
                p = p + h*n2; prev = n2
                if np.max(np.abs(p)) > box or np.hypot(*p) < rcore:
                    break
                pts.append(p.copy())
            lines.append(np.array(pts))
    return lines


def fig_director_lines(tagA):
    paper, full = load(tagA, 1), load(tagA, 4)
    if paper is None or full is None:
        return
    seeds = [(0.0, yy) for yy in np.arange(-9, 9.01, 1.5) if abs(yy) > 0.1] + \
            [(xx, 0.0) for xx in (-8, -5, -3, -1.5)]
    fig, ax = plt.subplots(figsize=(3.375, 3.2))
    fig.subplots_adjust(left=0.15, right=0.97, bottom=0.13, top=0.97)
    for d, col, ls, lw in ((paper, "black", "--", 0.7), (full, RED, "-", 0.8)):
        for ln in director_lines(d, seeds):
            if len(ln) > 1:
                ax.plot(ln[:, 0], ln[:, 1], ls, color=col, lw=lw)
    ax.plot([0], [0], "o", ms=3, color="k")
    ax.set_aspect("equal")
    ax.set_xlim(-9.5, 9.5); ax.set_ylim(-9.5, 9.5)
    ax.set_xlabel(r"$x$ (%s)" % UM, labelpad=1)
    ax.set_ylabel(r"$y$ (%s)" % UM)
    save(fig, "director_lines.pdf")


# ---------------------------------------------------------------- u against the box
def fig_u_box(box_tags):
    rows = []
    for L, tag in box_tags:
        r = load_runs(tag)
        if r is None:
            continue
        d = {tuple(np.round(m["w"], 6)): (m["u"][0] if m["ok"] else np.nan)
             for m in r["milestones"]}
        rows.append((L, d.get((0, 0, 0, 0), np.nan), d.get((1, 0, 0, 0), np.nan),
                     d.get((1, 1, 0, 0), np.nan), d.get((1, 1, 1, 1), np.nan)))
    if len(rows) < 2:
        return
    rows = np.array(sorted(rows))
    fig, ax = plt.subplots(figsize=(3.375, 2.5))
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.96)
    ax.axhline(0, color=GRAY, lw=0.6, zorder=0)
    for j, lab, mk, col, mfc in ((1, "paper", "o", "black", "white"),
                                 (2, r"$+w_T$", "s", ORANGE, "white"),
                                 (3, r"$+w_R$", "^", ORANGE, ORANGE),
                                 (4, "full", "o", BLUE, BLUE)):
        ax.plot(rows[:, 0], rows[:, j], "-", color=col, lw=0.8, zorder=1)
        ax.plot(rows[:, 0], rows[:, j], mk, ms=4, mfc=mfc, mec=col, color=col, zorder=2,
                label=lab)
    ax.set_xlabel(r"half-width $L$ of the box (%s)" % UM, labelpad=1)
    ax.set_ylabel(r"$u_x$ (%s/min)" % UM)
    ax.legend(loc="best", frameon=False, handletextpad=0.3, ncol=2)
    save(fig, "u_box.pdf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--A", default="A1")
    ap.add_argument("--B", default="B1")
    a = ap.parse_args()
    fig_u_path(a.A, a.B)
    fig_S_axes(a.A)
    fig_director_dev(a.A)
    fig_rho_axis(a.A)
    fig_forces_axis(a.A)
    fig_maps(a.A)
    fig_director_lines(a.A)
    fig_u_box([(10, "L10"), (15, "L15"), (20, a.A), (30, "L30"), (40, "L40")])


if __name__ == "__main__":
    main()
