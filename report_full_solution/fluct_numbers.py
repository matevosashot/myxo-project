#!/usr/bin/env python3
"""Numbers of the fluctuation section: data/fluct_summary.json (fluct_analysis.py)
and the check logs of ../lyapunov_full/logs/ -> numbers_fluct.tex (macros \\nf...).
Run: ../.venv/bin/python -I fluct_numbers.py [--N 180]
"""
import argparse
import json
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LOGS = os.path.join(HERE, "..", "lyapunov_full", "logs")
lines = []


def put(name, value, fmt="{}"):
    lines.append("\\newcommand{\\nf%s}{%s}" % (name, fmt.format(value)))


def sci(v, d=1):
    m, e = ("%.*e" % (d, v)).split("e")
    return r"%s\times10^{%d}" % (m, int(e))


def pick(summ, w, N=None, box=20):
    c = [r for r in summ["runs"].values() if r["w"] == w and r["maxReLambda"] < 0
         and r["fbox"] == box and (N is None or r["N"] == N)]
    return max(c, key=lambda r: r["N"]) if c else None


def logval(fname, pat):
    p = os.path.join(LOGS, fname)
    if not os.path.exists(p):
        return None
    m = re.search(pat, open(p).read())
    return float(m.group(1).replace("*^", "e")) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=None)
    N = ap.parse_args().N
    summ = json.load(open(os.path.join(HERE, "data", "fluct_summary.json")))
    if N is None:
        ok = [n for n in sorted({r["N"] for r in summ["runs"].values() if r["fbox"] == 20})
              if pick(summ, "0000", n) and pick(summ, "1111", n)]
        N = ok[-1]
    put("N", N)
    P, F = pick(summ, "0000", N), pick(summ, "1111", N)
    put("h", P["h"], "{:.3f}")
    put("lEff", P["lEff"], "{:.3f}")
    put("ped", sci(P["pedestal"], 2))
    for lbl, r in (("Paper", P), ("Full", F)):
        put("nu" + lbl, 100 * r["nu0"], "{:.2f}")
        put("x" + lbl, r["x0"], "{:.2f}")
        put("nuMin" + lbl, 100 * r["nuMin"], "{:.2f}")
        put("xMin" + lbl, r["xMin"], "{:.2f}")
        put("nuCross" + lbl, 100 * r["nuCrossMax"], "{:.2f}")
        put("R" + lbl, sci(r["R"], 1) if r["R"] >= 1e3 else "%.0f" % r["R"])
        put("nuAn" + lbl, 100 * r["nu0Analytic"], "{:.2f}")
        put("maxRe" + lbl, sci(r["maxReLambda"], 1))
        put("res" + lbl, sci(r["residual"], 0))
        put("sigS" + lbl, r["sigmaS_far"], "{:.2f}")
        put("sigT" + lbl, r["sigmaTheta_far_deg"], "{:.0f}")
        put("coreX" + lbl, r["coreStd"][0], "{:.1f}")
        put("coreY" + lbl, r["coreStd"][1], "{:.1f}")
        put("slow" + lbl, -r["slowest"][0][0], "{:.4f}")
        put("slowB" + lbl, -r["slowest"][1][0], "{:.4f}")
        put("trans" + lbl, min(r["translationOverlap"][:2]), "{:.2f}")
        put("srcsum" + lbl, sci(r["sourcesSumVsTotal"][0], 0))
        sm = r["src_max"]
        put("aniso" + lbl, 100 * sm["rhoAniso"], "{:.2f}")
        for k in ("cross", "QQrev", "QQrot"):
            put(k + lbl, 1e4 * sm[k], "{:.1f}")
        pm = r["path_max"]
        put("pathQ" + lbl, 1e4 * pm["nematic"], "{:.1f}")
        put("pathI" + lbl, 1e4 * pm["interference"], "{:.1f}")
        put("pathDirect" + lbl, 100 * (pm["direct"] - 1), "{:.2f}")
        rng = r["srcRange"]
        put("QmaxAbs" + lbl, 1e4 * max(abs(v) for k in ("cross", "QQrev", "QQrot")
                                        for v in rng[k]), "{:.0f}")
        put("pathQmax" + lbl, 1e4 * max(abs(v) for v in r["pathRange"]["nematic"]), "{:.0f}")
        sf = r["srcFar"]
        put("isoFar" + lbl, 100 * (sf["rhoIso"] - 1), "{:.2f}")
        put("crossFar" + lbl, 100 * sf["cross"], "{:.2f}")
        put("QQrevFar" + lbl, 100 * sf["QQrev"], "{:.2f}")
        put("QQrotFar" + lbl, 1e4 * sf["QQrot"], "{:.1f}")
        put("pathQFar" + lbl, 1e4 * r["pathFar"]["nematic"], "{:.1f}")
        put("shareRot" + lbl, 100 * r["nematicShareFar"]["QQrot"][1], "{:.0f}")
    # ladder on the finest grid where all four are stable
    for w, lbl in (("1000", "WT"), ("1100", "WTR")):
        r = pick(summ, w)
        if r:
            put("nu" + lbl, 100 * r["nu0"], "{:.2f}")
            put("x" + lbl, r["x0"], "{:.2f}")
            put("R" + lbl, sci(r["R"], 1))
            put("N" + lbl, r["N"])
    lad = pick(summ, "0000", pick(summ, "1000")["N"]) if pick(summ, "1000") else None
    if lad:
        put("nuPaperLadder", 100 * lad["nu0"], "{:.2f}")
    # reference: the paper's ssDir runs
    for n, rr in summ["ref"].items():
        put("nuRef" + ("A" if n == "150" else "B"), 100 * rr["nu0"], "{:.2f}")
        put("xRef" + ("A" if n == "150" else "B"), rr["x0"], "{:.2f}")
        put("RRef" + ("A" if n == "150" else "B"), "%.0f" % rr["R"])
        if "maxDiffVsP0" in rr:
            put("refDiff" + ("A" if n == "150" else "B"), sci(rr["maxDiffVsP0"], 0))
    # box scan
    put("hBox", 40.0 / 151, "{:.3f}")
    for box in (10, 15, 20, 30):
        for w, lbl in (("0000", "Paper"), ("1111", "Full")):
            c = [r for r in summ["runs"].values() if r["w"] == w and r["fbox"] == box
                 and r["maxReLambda"] < 0 and abs(r["h"] / (40.0 / 151) - 1) < 0.03]
            r = c[0] if c else None
            if r:
                tag = {10: "Ten", 15: "Fifteen", 20: "Twenty", 30: "Thirty"}[box]
                put("nu" + lbl + tag, 100 * r["nu0"], "{:.2f}")
                put("slow" + lbl + tag, -r["slowest"][0][0], "{:.4f}")
                put("core" + lbl + tag, float(np.hypot(*r["coreStd"])), "{:.1f}")
                put("sigT" + lbl + tag, r["sigmaTheta_far_deg"], "{:.0f}")
    # checks
    put("chkA", logval("check_small_N41_m1.log", r"T1  \|A - A_module\|/\|A\|\s*=\s*(\S+)") or 0, "{:.0g}")
    v = logval("check_small_N41_m1.log", r"T1  \|WNW - Q_module\|/\|Q\|\s*=\s*(\S+)")
    put("chkQ", sci(v, 0) if v else "--")
    for m in ("1", "4"):
        v = logval("check_small_N41_m%s.log" % m, r"T5  \|C_py - C_ref\|/max\|C_ref\| = (\S+)")
        put("chkC" + ("P" if m == "1" else "F"), sci(v, 0) if v else "--")
        for key, pat in (("Src", r"T3  QQrot.*?= (\S+)"), ("Vd", r"T4  var_d\s+rel.err = (\S+)"),
                         ("Cd", r"T4  cov_d\s+rel.err = (\S+)")):
            v = logval("check_multi_N41_m%s.log" % m, pat)
            put("chk" + key + ("P" if m == "1" else "F"), sci(v, 0) if v else "--")
    p = os.path.join(LOGS, "check_linearization.log")
    if os.path.exists(p):
        txt = open(p).read()
        for nI, lbl in ((41, "A"), (81, "B"), (161, "C")):
            m = re.search(r"Nint %d  w = \{1, 1, 1, 1\}.*?= \{([^}]*)\}" % nI, txt)
            if m:
                put("lin" + lbl, max(float(x) for x in m.group(1).split(",")), "{:.3f}")
    open(os.path.join(HERE, "numbers_fluct.tex"), "w").write(
        "% generated by fluct_numbers.py -- do not edit\n" + "\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
