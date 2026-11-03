#!/usr/bin/env python3
"""Figures of the fluctuation section of report.tex.

Data:   data/fluct_summary.json and data/fluct_<tag>.npz (fluct_analysis.py),
        data/fluct_ref_N<N>.json (export_fluct_ref.wls, the paper's ssDir runs).
Run:    ../.venv/bin/python -I make_fluct_figures.py [--N 180]
        --N picks the grid of the paper/full comparison (default: the finest grid
        on which both are stable).
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import (BLUE, ORANGE, RED, GREEN, GRAY, UM, DATA, panel_label, save)

LN = 0.5
Z_INF = 20.8
MODELS = [("0000", "paper", "k"), ("1000", r"$+w_T$", ORANGE),
          ("1100", r"$+w_T,w_R$", GREEN), ("1111", "full", BLUE)]
SRC = [("rhoIso", r"$\rho\rho$ isotropic $-1$", GRAY),
       ("rhoAniso", r"$\rho\rho$ anisotropic", BLUE),
       ("cross", r"$\rho$--$\mathbf Q$ cross", ORANGE),
       ("QQrev", r"$\mathbf{QQ}$ reversal", GREEN),
       ("QQrot", r"$\mathbf{QQ}$ rotational", RED)]


def load_npz(tag):
    return np.load(os.path.join(DATA, "fluct_%s.npz" % tag))


def axis(F, N):
    F = F.reshape(N, N)
    return 0.5 * (F[N // 2 - 1] + F[N // 2])


def cross(F, N):
    F = F.reshape(N, N)
    return 0.5 * (F[:, N // 2 - 1] + F[:, N // 2])


def pick(summ, w, N=None, box=20):
    c = [r for r in summ["runs"].values() if r["w"] == w and r["maxReLambda"] < 0
         and r["fbox"] == box and (N is None or r["N"] == N)]
    return max(c, key=lambda r: r["N"]) if c else None


def fig_axis(summ, N):
    """nu along and across the axis, paper vs full, Lyapunov and analytic"""
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.5))
    for w, lbl, col in (MODELS[0], MODELS[3]):
        r = pick(summ, w, N)
        if r is None:
            continue
        o = load_npz(r["tag"])
        g = o["g"]
        axs[0].plot(g, 100 * axis(o["nu"], r["N"]), color=col, lw=1.2, label=lbl)
        axs[0].plot(g, 100 * axis(o["nu_an"], r["N"]), color=col, lw=0.8, ls="--")
        axs[1].plot(g, 100 * cross(o["nu"], r["N"]), color=col, lw=1.2, label=lbl)
        axs[1].plot(g, 100 * cross(o["nu_an"], r["N"]), color=col, lw=0.8, ls="--")
    ref = os.path.join(DATA, "fluct_ref_N%d.json" % N)
    if os.path.exists(ref):
        rr = json.load(open(ref))
        ped = 12.5 / (2 * np.pi * 9e4 * 0.25 * LN**2)
        nu = np.array(rr["sigmaRho"]).ravel() / ped - 1
        axs[0].plot(rr["gridInt"], 100 * axis(nu, N), color="k", lw=0, marker=".", ms=2,
                    label=r"paper run (ssDir)")
    for ax, lab in zip(axs, ("x", "y")):
        ax.set_xlim(-7, 7)
        ax.axhline(0, color="0.8", lw=0.5)
        ax.set_xlabel(r"$%s$ (%s)" % (lab, UM))
    for ax in axs:
        ax.set_ylim(-4.5, 7)
    axs[0].set_ylabel(r"$\nu$ (\%)")
    axs[0].set_title(r"along the axis, $y=0$", fontsize=9)
    axs[1].set_title(r"across, $x=0$", fontsize=9)
    axs[0].legend(frameon=False, loc="upper right")
    panel_label(axs[0], "a"); panel_label(axs[1], "b")
    fig.tight_layout()
    save(fig, "fluct_axis.pdf")


def fig_maps(summ, N):
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.9))
    vals = []
    for w, lbl, col in (MODELS[0], MODELS[3]):
        r = pick(summ, w, N)
        if r is not None:
            o = load_npz(r["tag"])
            X, Y = np.meshgrid(o["g"], o["g"])
            vals.append(100 * o["nu"][(np.hypot(X, Y) < 5).ravel()])
    vmax = max(np.max(np.abs(v)) for v in vals) if vals else 1
    for ax, (w, lbl, col), s in zip(axs, (MODELS[0], MODELS[3]), "ab"):
        r = pick(summ, w, N)
        if r is None:
            continue
        o = load_npz(r["tag"])
        g = o["g"]
        Nn = r["N"]
        im = ax.pcolormesh(g, g, 100 * o["nu"].reshape(Nn, Nn), cmap="RdBu_r",
                           vmin=-vmax, vmax=vmax, shading="auto", rasterized=True)
        ax.set_xlim(-5, 5); ax.set_ylim(-5, 5); ax.set_aspect("equal")
        ax.set_title(lbl, fontsize=9)
        ax.set_xlabel(r"$x$ (%s)" % UM)
        ax.plot(r["x0"], 0, "k+", ms=5)
        panel_label(ax, s, x=-0.2)
    axs[0].set_ylabel(r"$y$ (%s)" % UM)
    cb = fig.colorbar(im, ax=axs, shrink=0.85, pad=0.02)
    cb.set_label(r"$\nu$ (\%)")
    save(fig, "fluct_maps.pdf")


def fig_sources(summ, N):
    """contributions to Var drho along the axis: by noise source and by pathway"""
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.7), sharey=True)
    for ax, (w, lbl, col), s in zip(axs, (MODELS[0], MODELS[3]), "ab"):
        r = pick(summ, w, N)
        if r is None:
            continue
        o = load_npz(r["tag"])
        g, Nn = o["g"], r["N"]
        A = lambda k: axis(o[k], Nn)
        ax.plot(g, A("src_rhoIso") - 1, color=GRAY, lw=1,
                label=r"reversal, $\rho\rho$ isotropic $-1$")
        ax.plot(g, A("src_rhoAniso"), color=BLUE, lw=1, label=r"reversal, $\rho\rho$ anisotropic")
        ax.plot(g, A("src_cross") + A("src_QQrev"), color=ORANGE, lw=1,
                label=r"reversal, $\rho$--$\mathbf Q$ and $\mathbf{QQ}$")
        ax.plot(g, A("src_QQrot"), color=RED, lw=1, label="rotational")
        ax.plot(g, A("var_Q"), color=GREEN, lw=1, ls="--", label=r"pathway via $\delta\mathbf Q$")
        ax.plot(g, A("nu"), color="k", lw=0.9, ls=(0, (1, 1.5)), label=r"total, $\nu$")
        ax.set_yscale("symlog", linthresh=1e-4)
        ax.set_xlim(-7, 7)
        ax.set_xlabel(r"$x$ (%s)" % UM)
        ax.set_title(lbl + r", $y=0$", fontsize=9)
        panel_label(ax, s)
    axs[0].set_ylabel(r"contribution to $\mathrm{Var}\,\delta\rho/\Sigma_\mathrm{ped}$")
    axs[0].legend(frameon=False, fontsize=6, loc="lower left")
    fig.tight_layout()
    save(fig, "fluct_sources.pdf")


def fig_ladder(summ):
    """nu0 and the enhancement along the model ladder and against the grid"""
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.5))
    Ns = sorted({r["N"] for r in summ["runs"].values()})
    for w, lbl, col in MODELS:
        pts = sorted((r["h"], 100 * r["nu0"]) for r in summ["runs"].values()
                     if r["w"] == w and r["maxReLambda"] < 0 and r["fbox"] == 20)
        if pts:
            h, v = zip(*pts)
            axs[0].plot(np.array(h)**2, v, "o-", color=col, ms=3, lw=0.8, label=lbl)
    for N, rr in summ["ref"].items():
        h = 40 / (int(N) + 1)
        axs[0].plot(h**2, 100 * rr["nu0"], "kx", ms=5)
    axs[0].set_xlim(left=0)
    axs[0].set_xlabel(r"$\Delta x^2$ (%s$^2$)" % UM)
    axs[0].set_ylabel(r"$\nu_0$ (\%)")
    axs[0].legend(frameon=False, fontsize=7)
    nu = np.linspace(0, 0.07, 200)
    axs[1].semilogy(100 * nu, np.exp(0.5 * Z_INF**2 * nu / (1 + nu)), color=GRAY, lw=0.8)
    for w, lbl, col in MODELS:
        r = pick(summ, w)
        if r is not None:
            axs[1].semilogy(100 * r["nu0"], r["R"], "o", color=col, ms=4)
            axs[1].annotate(lbl, (100 * r["nu0"], r["R"]), xytext=(4, -8),
                            textcoords="offset points", fontsize=7, color=col)
    axs[1].set_xlabel(r"$\nu_0$ (\%)")
    axs[1].set_ylabel(r"$\mathcal P(\mathbf r_\mathrm{max})/\mathcal P(\infty)$")
    panel_label(axs[0], "a"); panel_label(axs[1], "b")
    fig.tight_layout()
    save(fig, "fluct_ladder.pdf")


def circle_profile(g, nu, center, rs):
    from scipy.interpolate import RegularGridInterpolator
    N = len(g)
    P = np.exp(0.5 * Z_INF**2 * (1 - 1 / (1 + nu.reshape(N, N))))
    f = RegularGridInterpolator((g, g), P.T)
    th = np.linspace(0, 2 * np.pi, 361)[:-1]
    return np.array([np.mean(f(np.stack([center[0] + r * np.cos(th),
                                         center[1] + r * np.sin(th)], -1))) for r in rs])


def fig_layering(summ, N):
    """circle-averaged layering probability, eq:p-layer-radial, about r_max and the core"""
    fig, ax = plt.subplots(figsize=(3.3, 2.5))
    rs = np.linspace(0, 6, 61)
    for w, lbl, col in (MODELS[0], MODELS[3]):
        r = pick(summ, w, N)
        if r is None:
            continue
        o = load_npz(r["tag"])
        ax.semilogy(rs, circle_profile(o["g"], o["nu"], (r["x0"], 0.0), rs), color=col,
                    lw=1.2, label=lbl + r", about $\mathbf r_\mathrm{max}$")
        ax.semilogy(rs, circle_profile(o["g"], o["nu"], (0.0, 0.0), rs), color=col,
                    lw=0.8, ls="--", label=lbl + ", about the core")
    ax.set_xlabel(r"$r$ (%s)" % UM)
    ax.set_ylabel(r"$\mathcal P(r)/\mathcal P(\infty)$")
    ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout()
    save(fig, "fluct_layering.pdf")


def fig_nematic(summ, N):
    """order and director fluctuations, and the density-order correlation"""
    fig, axs = plt.subplots(1, 3, figsize=(6.8, 2.4))
    for w, lbl, col in (MODELS[0], MODELS[3]):
        r = pick(summ, w, N)
        if r is None:
            continue
        o = load_npz(r["tag"])
        g, Nn = o["g"], r["N"]
        far = np.abs(g) > 2.0           # inside, the translation mode dominates dQ
        for side in (g < -2.0, g > 2.0):
            axs[0].plot(g[side], np.sqrt(np.clip(axis(o["varS"], Nn), 0, None))[side],
                        color=col, lw=1.1, label=lbl if side[0] else None)
            axs[1].plot(g[side], np.degrees(np.sqrt(np.clip(axis(o["varT"], Nn), 0, None)))[side],
                        color=col, lw=1.1)
            axs[2].plot(g[side], axis(o["corrRS"], Nn)[side], color=col, lw=1.1)
    axs[0].set_ylabel(r"$\sigma_S$")
    axs[1].set_ylabel(r"$\sigma_\theta$ (deg)")
    axs[2].set_ylabel(r"corr$(\delta\rho,\delta S)$")
    for ax, s in zip(axs, "abc"):
        ax.set_xlim(-10, 10)
        ax.set_xlabel(r"$x$ (%s)" % UM)
        panel_label(ax, s, x=-0.3)
    axs[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, "fluct_nematic.pdf")


def fig_corr(summ, N):
    """equal-time density correlation coefficient with the point of maximum variance"""
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.9))
    for ax, (w, lbl, col), s in zip(axs, (MODELS[0], MODELS[3]), "ab"):
        r = pick(summ, w, N)
        if r is None:
            continue
        o = load_npz(r["tag"])
        g, Nn = o["g"], r["N"]
        k = -2                                   # the row at argmax Var drho (lyap_multi.py)
        x0, y0 = o["row_xy"][k]
        k0 = int(np.argmin(np.abs(g - y0))) * Nn + int(np.argmin(np.abs(g - x0)))
        row, var = o["rows_rho"][k], o["var_tot"]
        c = row / np.sqrt(var[k0] * var)
        from matplotlib.colors import TwoSlopeNorm
        im = ax.pcolormesh(g, g, c.reshape(Nn, Nn), cmap="RdBu_r",
                           norm=TwoSlopeNorm(vmin=-0.1, vcenter=0, vmax=1),
                           shading="auto", rasterized=True)
        ax.plot(*o["row_xy"][k], "k+", ms=5)
        ax.set_xlim(-4, 4); ax.set_ylim(-4, 4); ax.set_aspect("equal")
        ax.set_title(lbl, fontsize=9)
        ax.set_xlabel(r"$x$ (%s)" % UM)
        panel_label(ax, s, x=-0.2)
    axs[0].set_ylabel(r"$y$ (%s)" % UM)
    cb = fig.colorbar(im, ax=axs, shrink=0.85, pad=0.02)
    cb.set_label(r"$\langle\delta\rho(\mathbf r_\mathrm{max})\delta\rho(\mathbf r)\rangle$, normalized")
    save(fig, "fluct_corr.pdf")


def fig_box(summ):
    """box size: slowest relaxation rate (defect translation), nu0, core spread"""
    fig, axs = plt.subplots(1, 3, figsize=(6.8, 2.3))
    href = 40.0 / 151        # the N_int = 150 spacing, run for every box
    for w, lbl, col in (MODELS[0], MODELS[3]):
        pts = []
        for box in (10, 15, 20, 30):
            c = [r for r in summ["runs"].values() if r["w"] == w and r["fbox"] == box
                 and r["maxReLambda"] < 0 and abs(r["h"] / href - 1) < 0.03]
            r = c[0] if c else None
            if r:
                pts.append((box, -r["slowest"][0][0], 100 * r["nu0"],
                            np.hypot(*r["coreStd"])))
        if not pts:
            continue
        L, rate, nu0, core = map(np.array, zip(*pts))
        axs[0].loglog(L, rate, "o-", color=col, ms=3, lw=0.8, label=lbl)
        axs[1].plot(L, nu0, "o-", color=col, ms=3, lw=0.8)
        axs[2].plot(L, core, "o-", color=col, ms=3, lw=0.8)
    axs[0].set_ylabel(r"slowest rate (min$^{-1}$)")
    axs[1].set_ylabel(r"$\nu_0$ (\%)")
    axs[1].set_ylim(bottom=0)
    axs[2].set_ylabel(r"core spread (%s)" % UM)
    for ax, s in zip(axs, "abc"):
        ax.set_xlabel(r"box half-width $L$ (%s)" % UM)
        panel_label(ax, s, x=-0.3)
    axs[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, "fluct_box.pdf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=None)
    a = ap.parse_args()
    summ = json.load(open(os.path.join(DATA, "fluct_summary.json")))
    N = a.N
    if N is None:
        ok = [n for n in sorted({r["N"] for r in summ["runs"].values()})
              if pick(summ, "0000", n) and pick(summ, "1111", n)]
        N = ok[-1] if ok else max(r["N"] for r in summ["runs"].values())
    print("comparison grid N =", N)
    fig_axis(summ, N)
    fig_maps(summ, N)
    fig_sources(summ, N)
    fig_ladder(summ)
    fig_layering(summ, N)
    fig_nematic(summ, N)
    fig_corr(summ, N)
    fig_box(summ)


if __name__ == "__main__":
    main()
