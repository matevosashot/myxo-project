#!/usr/bin/env python3
"""fluct_analysis.py -- numbers of the fluctuation section from the runs of
../lyapunov_full/run_fluct.wls (../data/fluctuations-full/<tag>/) and the paper's
ssDir runs (data/fluct_ref_N<N>.json, from export_fluct_ref.wls).

Writes data/fluct_summary.json: per run the axis cuts, nu0 and its position, the
enhancement of eq:R-value, the source and pathway splits at the extrema, the
nematic fluctuations, the core-position covariance, the slowest relaxation
rates, and the analytic eq:total-gauss evaluated on the same steady state.
Field maps for the figures go to data/fluct_<tag>.npz.

Conventions: Var drho is divided by the pedestal zeta/(2 pi B rho0 ln^2); nu =
Var/pedestal - 1.  Even Nint puts no node on y = 0: the axis is the mean of the
two rows y = +-h/2 (the variance is even in y).
"""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "..", "data", "fluctuations-full")
Z_INF = 20.8          # manuscript eq:R-value (B = 9e4, ln = 0.5)
LN = 0.5
SRC = ["rhoIso", "rhoAniso", "cross", "QQrev", "QQrot"]
MODEL = {"0000": "paper", "1000": "+wT", "1100": "+wT,wR", "1111": "full"}


def enhancement(nu):
    return float(np.exp(0.5 * Z_INF**2 * nu / (1.0 + nu)))


def axis(field, N):
    """mean of the rows j = N/2-1, N/2 (y = -h/2, +h/2) of a [j, i] map"""
    F = field.reshape(N, N)
    return 0.5 * (F[N // 2 - 1] + F[N // 2]) if N % 2 == 0 else F[N // 2]


def cross_axis(field, N):
    F = field.reshape(N, N)
    return 0.5 * (F[:, N // 2 - 1] + F[:, N // 2]) if N % 2 == 0 else F[:, N // 2]


def extrema(x, nu, rmax=5.0):
    m = np.abs(x) < rmax
    i = np.flatnonzero(m)
    return i[np.argmax(nu[m])], i[np.argmin(nu[m])]


def analytic_nu(Q1, Q2, g, ln=LN):
    """(ln^2/4) (kappa * d_i d_j Q_ij)(r), eq:total-gauss, on the node grid (FFT)"""
    N = len(g)
    h = g[1] - g[0]
    A = Q1.reshape(N, N)
    B = Q2.reshape(N, N)
    P = lambda F: np.pad(F, 1)                      # Q = 0 beyond the wall nodes
    A, B = P(A), P(B)
    dxx = (A[1:-1, 2:] - 2 * A[1:-1, 1:-1] + A[1:-1, :-2]) / h**2
    dyy = (A[2:, 1:-1] - 2 * A[1:-1, 1:-1] + A[:-2, 1:-1]) / h**2
    dxy = (B[2:, 2:] - B[2:, :-2] - B[:-2, 2:] + B[:-2, :-2]) / (4 * h**2)
    f = dxx - dyy + 2 * dxy
    k = np.arange(-(N - 1), N) * h
    X, Y = np.meshgrid(k, k)
    w2 = X**2 + Y**2
    with np.errstate(invalid="ignore", divide="ignore"):
        kap = (1 - np.exp(-2 * w2 / ln**2)) / (2 * np.pi * w2)
    kap[w2 == 0] = 1.0 / (np.pi * ln**2)
    M = 2 * N - 1 + N
    conv = np.fft.irfft2(np.fft.rfft2(f, (M, M)) * np.fft.rfft2(kap, (M, M)), (M, M))
    conv = conv[N - 1:2 * N - 1, N - 1:2 * N - 1] * h**2
    return (ln**2 / 4) * conv.ravel()


def nematic(Q1, Q2, d):
    """sigma_S^2, sigma_theta^2 and Cov(drho, dS) from the six block diagonals"""
    S = np.hypot(Q1, Q2)
    c11, c22, c12, cr1, cr2 = d[1], d[2], d[5], d[3], d[4]
    varS = (Q1**2 * c11 + Q2**2 * c22 + 2 * Q1 * Q2 * c12) / S**2
    varT = (Q1**2 * c22 + Q2**2 * c11 - 2 * Q1 * Q2 * c12) / (4 * S**4)
    covRS = (Q1 * cr1 + Q2 * cr2) / S
    return S, varS, varT, covRS


def run_summary(d):
    o = np.load(os.path.join(d, "out.npz"))
    st = json.load(open(os.path.join(d, "steady.json")))
    res = json.load(open(os.path.join(d, "result.json")))
    N, m = st["nInt"], st["nInt"] ** 2
    g = np.array(st["gridInt"])
    h = st["h"]
    ped = st["pedestal"]
    w = "".join(str(int(round(v))) for v in st["w"])
    F = {k: np.array(v) for k, v in st["fields"].items()}
    tot = o["total"]
    nu = tot[0] / ped - 1
    nu_ax, nu_cr = axis(nu, N), cross_axis(nu, N)
    imx, imn = extrema(g, nu_ax)
    jmx = np.argmax(np.abs(nu_cr) * (np.abs(g) < 5))
    out = {"tag": os.path.basename(d), "N": N, "fbox": st["fbox"], "h": h, "w": w, "model": MODEL.get(w, w),
           "u": st["u"][0], "pedestal": ped, "lEff": st["lEff"],
           "maxReLambda": res["max_re_lambda"], "residual": res["residual_refined"],
           "sourcesSumVsTotal": [float(v) for v in
               np.abs(sum(o["src_" + k] for k in SRC) - o["total"]).max(axis=1) /
               np.abs(o["total"]).max(axis=1)],          # per block rr,11,22,r1,r2,12
           "minEigLocal": st["minEigLocal"], "peclet": st["peclet"],
           "timings": res["timings"], "peakRSS": res.get("peak_rss_gib")}
    out["nu0"] = float(nu_ax[imx]); out["x0"] = float(g[imx])
    out["nuMin"] = float(nu_ax[imn]); out["xMin"] = float(g[imn])
    out["nuCrossMax"] = float(nu_cr[jmx]); out["yCrossMax"] = float(g[jmx])
    out["R"] = enhancement(out["nu0"])
    out["nuFar"] = float(np.mean(nu.reshape(N, N)[np.ix_(np.abs(g) < 12, np.abs(g) < 12)]
                                 [np.hypot(*np.meshgrid(g[np.abs(g) < 12],
                                                        g[np.abs(g) < 12])) > 8]))
    # sources and pathways at the axis extrema (in units of the pedestal)
    srcAx = {s: axis(o["src_" + s][0], N) / ped for s in SRC}
    var_d = o["var_d"] / ped
    cov = o["cov_d"][0] / ped
    var_Q = tot[0] / ped + var_d - 2 * cov
    inter = 2 * (cov - var_d)
    path = {"direct": axis(var_d, N), "nematic": axis(var_Q, N),
            "interference": axis(inter, N)}
    for lbl, i in (("max", imx), ("min", imn)):
        out["src_" + lbl] = {s: float(srcAx[s][i]) for s in SRC}
        out["path_" + lbl] = {k: float(v[i]) for k, v in path.items()}
    rr = np.hypot(*np.meshgrid(g, g)).ravel()
    farm = (rr > 6) & (rr < 12)
    out["srcFar"] = {s: float(np.median(o["src_" + s][0][farm] / ped)) for s in SRC}
    out["pathFar"] = {"nematic": float(np.median(var_Q[farm])),
                      "interference": float(np.median(inter[farm]))}
    out["srcRange"] = {s: [float((o["src_" + s][0] / ped).min()),
                           float((o["src_" + s][0] / ped).max())] for s in SRC}
    out["pathRange"] = {"nematic": [float(var_Q.min()), float(var_Q.max())],
                        "interference": [float(inter.min()), float(inter.max())],
                        "directMinus1": [float(var_d.min() - 1), float(var_d.max() - 1)]}
    # analytic eq:total-gauss on this steady state
    nu_an = analytic_nu(F["Q1"], F["Q2"], g)
    nan_ax = axis(nu_an, N)
    ia, _ = extrema(g, nan_ax)
    out["nu0Analytic"] = float(nan_ax[ia]); out["x0Analytic"] = float(g[ia])
    # nematic fluctuations
    S, varS, varT, covRS = nematic(F["Q1"], F["Q2"], tot)
    far = np.hypot(*np.meshgrid(g, g)).ravel()
    ring = (far > 6) & (far < 10)
    out["sigmaS_far"] = float(np.sqrt(np.mean(varS[ring])))
    out["sigmaTheta_far_deg"] = float(np.degrees(np.sqrt(np.mean(varT[ring]))))
    corr = covRS / np.sqrt(tot[0] * varS)
    out["corrRhoS_range"] = [float(np.nanmin(corr[far > 1])), float(np.nanmax(corr[far > 1]))]
    srcS = {}
    for s in SRC:
        _, vS, vT, _ = nematic(F["Q1"], F["Q2"], o["src_" + s])
        srcS[s] = [float(np.mean(vS[ring]) / np.mean(varS[ring])),
                   float(np.mean(vT[ring]) / np.mean(varT[ring]))]
    out["nematicShareFar"] = srcS              # [share of var S, share of var theta]
    # core position: q(r) ~ J r + dq(0), dq(0) = mean of the four nodes around 0
    ri = o["row_index"]; rows = o["rows"]
    core = ri[:4]
    Sq = np.zeros((2, 2))
    for p in range(4):
        for a in (1, 2):
            for b in (1, 2):
                Sq[a - 1, b - 1] += sum(rows[p][a][b * m + core[pp]] for pp in range(4)) / 16
    Q1g, Q2g = F["Q1"].reshape(N, N), F["Q2"].reshape(N, N)
    j0, i0 = N // 2 - 1, N // 2 - 1
    dx = lambda G: 0.5 * ((G[j0, i0 + 1] - G[j0, i0]) + (G[j0 + 1, i0 + 1] - G[j0 + 1, i0])) / h
    dy = lambda G: 0.5 * ((G[j0 + 1, i0] - G[j0, i0]) + (G[j0 + 1, i0 + 1] - G[j0, i0 + 1])) / h
    J = np.array([[dx(Q1g), dy(Q1g)], [dx(Q2g), dy(Q2g)]])
    Ji = np.linalg.inv(J)
    Sr = Ji @ Sq @ Ji.T
    out["coreStd"] = [float(np.sqrt(max(Sr[0, 0], 0))), float(np.sqrt(max(Sr[1, 1], 0)))]
    # slowest relaxation
    mv = o["modes_val"]
    out["slowest"] = [[float(z.real), float(z.imag)] for z in mv[:6]]
    # are the slowest modes translations?  overlap of their dQ part with the span of
    # d_x Q_ss and d_y Q_ss (central differences on the node grid)
    gx = [np.gradient(G, h, axis=1).ravel() for G in (Q1g, Q2g)]
    gy = [np.gradient(G, h, axis=0).ravel() for G in (Q1g, Q2g)]
    T = np.stack([np.concatenate(gx), np.concatenate(gy)], 1)
    Tq, _ = np.linalg.qr(T)
    ov = []
    for k in range(min(4, o["modes_vec"].shape[1])):
        v = o["modes_vec"][m:, k].astype(complex)
        ov.append(float(np.linalg.norm(Tq.T @ v) / np.linalg.norm(v)))
    out["translationOverlap"] = ov
    # maps for the figures
    np.savez(os.path.join(HERE, "data", "fluct_%s.npz" % out["tag"]),
             g=g, nu=nu, nu_an=nu_an, var_d=var_d, var_Q=var_Q, inter=inter,
             S=S, varS=varS, varT=varT, corrRS=corr,
             **{"src_" + s: o["src_" + s][0] / ped for s in SRC},
             rows_rho=np.array([rows[k][0][:m] for k in range(len(ri))]),
             row_xy=o["row_xy"], var_tot=tot[0],
             modes_val=mv, modes_vec=o["modes_vec"][:, :4])
    return out


def ref_summary(N):
    f = os.path.join(HERE, "data", "fluct_ref_N%d.json" % N)
    if not os.path.exists(f):
        return None
    r = json.load(open(f))
    g = np.array(r["gridInt"])
    ped = 12.5 / (2 * np.pi * 9e4 * 0.25 * LN**2)
    nu = np.array(r["sigmaRho"]).ravel() / ped - 1
    nu_ax = axis(nu, N)
    imx, imn = extrema(g, nu_ax)
    return {"N": N, "nu0": float(nu_ax[imx]), "x0": float(g[imx]),
            "nuMin": float(nu_ax[imn]), "xMin": float(g[imn]),
            "R": enhancement(float(nu_ax[imx])), "maxReLambda": r["maxReLambda"],
            "nu": nu}


def main():
    runs = sorted(glob.glob(os.path.join(RUNS, "box*_N*")))
    summ = {"runs": {}, "ref": {}}
    for d in runs:
        if not os.path.exists(os.path.join(d, "out.npz")):
            continue
        s = run_summary(d)
        summ["runs"][s["tag"]] = s
        print("%-40s N=%3d %-7s maxRe=%+.2e nu0=%+.4f at x=%+.2f  nuMin=%+.4f at %+.2f  "
              "R=%.3g  analytic %.4f" % (s["tag"], s["N"], s["model"], s["maxReLambda"],
                                         s["nu0"], s["x0"], s["nuMin"], s["xMin"], s["R"],
                                         s["nu0Analytic"]))
    for N in (150, 180):
        r = ref_summary(N)
        if r is None:
            continue
        nu = r.pop("nu")
        summ["ref"][str(N)] = r
        mine = [s for s in summ["runs"].values() if s["N"] == N and s["w"] == "0000"]
        if mine:
            o = np.load(os.path.join(HERE, "data", "fluct_%s.npz" % mine[0]["tag"]))
            r["maxDiffVsP0"] = float(np.max(np.abs(o["nu"] - nu)))
        print("ref N%d: nu0=%.4f at %.2f  R=%.3g" % (N, r["nu0"], r["x0"], r["R"]))
    json.dump(summ, open(os.path.join(HERE, "data", "fluct_summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
