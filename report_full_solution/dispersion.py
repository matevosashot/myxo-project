#!/usr/bin/env python3
"""Linear stability of the uniform nematic state of the full model (continuum).

Independent of the Mathematica solver.  Base state: rho = 1, q = S (director
along x), h = h0 = -(a + 2 b S^2) S, C0 = g2 S^2, with S from the uniform balance
    (4/xir) h0 (1 - wR g2 S^2) = 4 Lambda S.
Perturbations ~ exp(i k.r), k = k (cos al, sin al).  The H equation is algebraic
in Fourier space, dh1 = A1 dq1, dh2 = A2 dq2 with
    A1 = -(a + 6 b S^2 + K k^2),   A2 = -(a + 2 b S^2 + K k^2),
so the growth rates are the eigenvalues of a 3x3 matrix acting on (drho, dq1, dq2)
(equations of comoving_full_core.m with u = 0, written in components).

Large-k limit with the density slaved (B k^2 >> K k^4): the amplitude mode dq1
grows at rate (wH K/xi0)(2 S^2 - 1 - wE g2 S^2) k^4, i.e. the grad-H pair is
ill-posed with the lowest-order mobility (wE = 0) for S > 1/sqrt(2), and
well-posed with the exact one (wE = 1), whose coefficient is -2 Var(cos 2 phi).

Writes data/dispersion.json and figures/dispersion.pdf.
Run:  ../.venv/bin/python -I dispersion.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from scipy.special import iv

HERE = os.path.dirname(os.path.abspath(__file__))

# bare parameters of the paper's steady state (meta of the B = 9e4 N180 run)
P = dict(B=9.0e4, a=-20.7, b=11.574074074074078, K=9.375, zeta=12.5, xi0=5.0,
         xir=15.0, Lam=0.13)


def kappa_of(S):
    return brentq(lambda k: iv(1, k)/iv(0, k) - S, 1e-12, 700.0)


def g2_and_slope(S):
    """g2(S) and dg2/ds2 (s2 = S^2) of the von Mises distribution."""
    k = kappa_of(S)
    m1 = iv(1, k)/iv(0, k); m2 = iv(2, k)/iv(0, k)
    dm1 = 1 - m1/k - m1**2
    dm2 = m1 - 2*m2/k - m1*m2
    g2 = m2/m1**2
    dg2 = (dm2/m1**2 - 2*m2*dm1/m1**3)/(2*m1*dm1)
    return g2, dg2


def base_state(w):
    wT, wR, wH, wE = w
    a, b, xir, Lam = P["a"], P["b"], P["xir"], P["Lam"]
    f = lambda S: -(a + 2*b*S*S)*(1 - wR*g2_and_slope(S)[0]*S*S) - Lam*xir
    S = brentq(f, 0.3, 0.995)
    g2, dg2 = g2_and_slope(S)
    h0 = -(a + 2*b*S*S)*S
    C0 = g2*S*S
    Cq = S*(dg2*S*S + 2*g2)          # dC/dq at q = S
    Cqb = dg2*S**3                   # dC/dqbar
    return dict(S=S, g2=g2, dg2=dg2, h0=h0, C0=C0, Cq=Cq, Cqb=Cqb)


def matrix(k, al, w, bs):
    wT, wR, wH, wE = w
    B, a, b, K, ze, x0, xr, Lam = (P[s] for s in ("B", "a", "b", "K", "zeta", "xi0", "xir", "Lam"))
    S, h0, C0, Cq, Cqb = (bs[s] for s in ("S", "h0", "C0", "Cq", "Cqb"))
    kx, ky = k*np.cos(al), k*np.sin(al)
    k2 = k*k; c2 = kx*kx - ky*ky; s2 = 2*kx*ky
    A1 = -(a + 6*b*S*S + K*k2)
    A2 = -(a + 2*b*S*S + K*k2)
    lam, rot = ze/x0, 4/xr
    M = np.zeros((3, 3))
    # rho row: (B/xi0) Lap rho + (zeta/xi0) Re dbar^2 q - (2 wH/xi0) Re(q Lap hbar + ...)
    M[0, 0] = -(B/x0)*k2
    M[0, 1] = -lam*c2 + (2*wH/x0)*S*k2*A1
    M[0, 2] = -lam*s2
    # q1 row
    M[1, 0] = -(B/x0)*S*k2
    M[1, 1] = ((wH/x0)*k2*(1 + wE*C0)*A1 + rot*((1 - wR*C0)*A1 - wR*h0*(Cq + Cqb))
               - 4*Lam - lam*k2 - lam*(wT/2)*c2*(Cq + Cqb))
    M[1, 2] = -lam*(wT/2)*s2*(Cq - Cqb)
    # q2 row
    M[2, 0] = 0.0
    M[2, 1] = lam*(wT/2)*s2*(Cq + Cqb)
    M[2, 2] = ((wH/x0)*k2*(1 - wE*C0)*A2 + rot*((1 + wR*C0)*A2 - wR*h0*(Cq - Cqb))
               - 4*Lam - lam*k2 - lam*(wT/2)*c2*(Cq - Cqb))
    return M


def max_growth(k, w, bs, nal=73):
    best = -np.inf
    for al in np.linspace(0, np.pi, nal):
        ev = np.linalg.eigvals(matrix(k, al, w, bs))
        best = max(best, ev.real.max())
    return best


CASES = [((0, 0, 0, 0), "paper", "black", "--"),
         ((1, 1, 0, 0), r"$w_T,w_R$", "#d9600b", "-."),
         ((1, 1, 1, 0), r"$w_T,w_R,w_H$ ($w_E=0$)", "#b8322a", "-"),
         ((1, 1, 1, 1), "full", "#3f6fb0", "-")]


def main():
    ks = np.logspace(-2, 3, 400)
    out = {"params": P, "k": ks.tolist(), "cases": []}
    for w, name, col, ls in CASES:
        bs = base_state(w)
        g = np.array([max_growth(k, w, bs) for k in ks])
        # rotation (Goldstone) mode at k -> 0 must be marginal
        ev0 = np.linalg.eigvals(matrix(1e-9, 0.0, w, bs))
        amp4 = (w[2]*P["K"]/P["xi0"])*(2*bs["S"]**2 - 1 - w[3]*bs["g2"]*bs["S"]**2)
        out["cases"].append({"w": list(w), "name": name, "S": bs["S"], "g2": bs["g2"],
                             "maxGrowth": g.tolist(), "eig_k0": sorted(ev0.real.tolist()),
                             "amplitude_k4_coefficient": amp4})
        unst = ks[g > 1e-9]
        print(f"{name:28s} S = {bs['S']:.4f}  g2 = {bs['g2']:.4f}  "
              f"k^4 coeff (amplitude) = {amp4:+.4f}  eig(k->0) = "
              f"{np.round(sorted(ev0.real), 6)}  "
              f"unstable k: {(unst.min(), unst.max()) if len(unst) else 'none'}")
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    with open(os.path.join(HERE, "data", "dispersion.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    plt.rcParams.update({
        "text.usetex": True, "font.family": "serif", "font.size": 9,
        "axes.labelsize": 10, "legend.fontsize": 8,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": False, "ytick.right": False, "axes.linewidth": 0.7,
        "text.latex.preamble": r"\usepackage{upgreek}\usepackage{amsmath}"})
    fig, ax = plt.subplots(figsize=(3.375, 2.6))
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.16, top=0.96)
    for c, (w, name, col, ls) in zip(out["cases"], CASES):
        g = np.array(c["maxGrowth"])
        ax.plot(ks, np.sign(g)*np.log10(1 + np.abs(g)), ls, color=col, lw=1.1, label=name)
    ax.axhline(0, color="0.45", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlim(ks[0], ks[-1])
    ax.set_xlabel(r"$k$ ($\upmu$m$^{-1}$)", labelpad=1)
    ax.set_ylabel(r"$\mathrm{sgn}(\lambda)\log_{10}(1+|\lambda|)$, $\lambda$ in min$^{-1}$")
    ax.axvline(np.pi/0.0749, color="0.6", lw=0.6, ls=":")
    ax.text(np.pi/0.0749*1.08, 1.0, r"grid $\pi/h_\mathrm{min}$", rotation=90, fontsize=7,
            color="0.45", va="bottom")
    ax.legend(loc="lower left", frameon=False, handlelength=1.8)
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", "dispersion.pdf"))
    print("wrote figures/dispersion.pdf, data/dispersion.json")


if __name__ == "__main__":
    main()
