#!/usr/bin/env python3
"""Numbers quoted in report.tex, collected from the exported runs.

Reads data/<tag>_runs.json and data/<tag>_m<k>.json (export_fields.wls) and
data/check_jacobian.log, and writes
    numbers.tex          LaTeX macros \\num...{} used by report.tex
    data/numbers.json    the same numbers, machine readable
Every number in the report text comes from here or from a figure script.
Run:  ../.venv/bin/python -I collect_numbers.py
"""
import json
import os
import re

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def runs(tag):
    p = os.path.join(DATA, f"{tag}_runs.json")
    return json.load(open(p)) if os.path.exists(p) else None


def mile(tag, k):
    p = os.path.join(DATA, f"{tag}_m{k}.json")
    return json.load(open(p)) if os.path.exists(p) else None


def find_w(run, w):
    for m in run["milestones"]:
        if np.allclose(m["w"], w):
            return m
    return None


def pade(r, S0=0.9, ld=1.0):
    x = r/ld
    return S0*x*np.sqrt((0.34 + 0.07*x**2)/(1 + 0.41*x**2 + 0.07*x**4))


def r_half(r, S, Sref):
    """first radius where S reaches Sref/2 (linear interpolation)"""
    r, S = np.asarray(r), np.asarray(S)
    i = np.argmax(S >= 0.5*Sref)
    return r[i-1] + (0.5*Sref - S[i-1])*(r[i] - r[i-1])/(S[i] - S[i-1])


def axis_halves(m):
    x = np.array(m["axisX"]["x"]); S = np.array(m["axisX"]["S"])
    y = np.array(m["axisY"]["y"]); Sy = np.array(m["axisY"]["S"])
    Sw = m["diag"]["wallS"]["+y"]
    tail = r_half(x[x >= 0], S[x >= 0], Sw)
    head = r_half(-x[x <= 0][::-1], S[x <= 0][::-1], Sw)
    side = r_half(y[y >= 0], Sy[y >= 0], Sw)
    return dict(tail=tail, head=head, side=side, Swall=Sw)


def max_dev(m, r):
    for c in m["circles"]:
        if abs(c["r"] - r) < 1e-9:
            return float(np.degrees(np.max(np.abs(c["dev"]))))
    return None


macros = {}
num = {}


def sci(v, digits=1):
    """LaTeX scientific notation, e.g. 2.2\\times10^{-5}"""
    if v == 0:
        return "0"
    e = int(np.floor(np.log10(abs(v))))
    m = v/10**e
    if round(abs(m), digits) >= 10:
        m /= 10; e += 1
    ms = f"{m:.{digits}f}" if digits > 0 else f"{m:.0f}"
    return ("" if ms in ("1", "1.0") else ms + "\\times") + f"10^{{{e}}}"


def put(name, value, fmt):
    if value is None:
        return
    num[name] = value
    macros[name] = sci(value, int(fmt[3:])) if fmt.startswith("sci") else fmt.format(value)


W = {"P0": [0, 0, 0, 0], "P1": [1, 0, 0, 0], "P2": [1, 1, 0, 0], "Full": [1, 1, 1, 1],
     "BwR": [0, 1, 0, 0], "BwRH": [0, 1, 1, 1], "Hlow": [1, 1, 1, 0]}
LETTERS = {"A1": "AOne", "B1": "BOne", "C1": "COne", "A15": "AOneHalf", "A2": "ATwo",
           "Z1": "ZOne", "Z15": "ZOneHalf", "Z2": "ZTwo", "H1": "HOne",
           "L10": "LTen", "L15": "LFifteen", "L30": "LThirty", "L40": "LForty",
           "Z1c": "ZOneC", "F15": "FOneHalf", "F2": "FTwo", "FZ2": "FZTwo",
           "G2": "GTwo"}

for tag, L in LETTERS.items():
    r = runs(tag)
    if r is None:
        continue
    for wn, w in W.items():
        m = find_w(r, w)
        if m is None:
            continue
        put(f"u{L}{wn}", m["u"][0], "{:+.4f}")
        put(f"uFine{L}{wn}", m["u"][0], "{:+.7f}")
        put(f"ok{L}{wn}", int(bool(m["ok"])), "{:d}")
        put(f"uSci{L}{wn}", m["u"][0], "sci1")
        put(f"uy{L}{wn}", abs(m["u"][1]), "sci0")
        put(f"res{L}{wn}", m["residual"], "sci0")
        put(f"lf{L}{wn}", m["lineForce"], "sci0")
    put(f"time{L}", r["time"]/60.0, "{:.0f}")
    put(f"nx{L}", r["nx"], "{:d}")
    put(f"hmin{L}", r["hMin"], "{:.3f}")

# path A at refine 1: structure at every milestone (paper, +wT, +wR, full)
pa, fu = mile("A1", 1), mile("A1", 4)
if pa and fu:
    for name, m in (("Paper", pa), ("PlusT", mile("A1", 2)), ("PlusR", mile("A1", 3)),
                    ("Full", fu)):
        if m is None:
            continue
        h = axis_halves(m)
        put(f"rhalfTail{name}", h["tail"], "{:.2f}")
        put(f"rhalfHead{name}", h["head"], "{:.2f}")
        put(f"rhalfSide{name}", h["side"], "{:.2f}")
        put(f"Swall{name}", h["Swall"], "{:.3f}")
        put(f"SwallPx{name}", m["diag"]["wallS"]["+x"], "{:.3f}")
        put(f"SwallMx{name}", m["diag"]["wallS"]["-x"], "{:.3f}")
        lo, hi = m["diag"]["rhoRange"]
        put(f"rhoMin{name}", lo, "{:.2f}")
        put(f"rhoMax{name}", hi, "{:.2f}")
        for rr in (0.5, 1.0, 2.0, 5.0):
            put(f"dev{name}R{str(rr).replace('.', 'p')}", max_dev(m, rr), "{:.1f}")
        x = np.array(m["axisX"]["x"])
        fa = np.abs(np.array(m["axisX"]["fActive"])); fe = np.abs(np.array(m["axisX"]["fElastic"]))
        core = np.abs(x) < 3
        put(f"fActiveMax{name}", fa[core].max(), "{:.1f}")
        if m["w"][2] != 0:          # the elastic force acts only with wH
            put(f"fElasticMax{name}", fe[core].max(), "{:.1f}")
        i0 = int(np.argmin(np.abs(x)))
        put(f"fNetCore{name}", m["axisX"]["fActive"][i0] + m["axisX"]["fPressure"][i0]
            + (m["axisX"]["fElastic"][i0] if m["w"][2] != 0 else 0.0), "{:.1f}")
    put("SInfUniform", fu["diag"]["SInfUniform"], "{:.3f}")
    rp = brentq(lambda r: pade(r) - 0.45, 0.01, 10)
    put("rhalfPade", rp, "{:.2f}")

# grid convergence of the full u: Richardson extrapolation from refinement 1.5 and 2 (h^4)
f15, f2, a1 = runs("F15"), runs("F2"), runs("A1")
if f15 and f2 and a1:
    u15, u2 = find_w(f15, W["Full"])["u"][0], find_w(f2, W["Full"])["u"][0]
    u1 = find_w(a1, W["Full"])["u"][0]
    uinf = u2 + (u2 - u15)/((2/1.5)**4 - 1)
    put("uExtrap", uinf, "{:+.5f}")
    put("uErrOne", abs(u1 - uinf), "sci1")
    put("uErrOneRel", 100*abs(u1 - uinf)/abs(uinf), "{:.2f}")

# the wE = 0 continuation: furthest wH reached (converged trace point at zeta target)
h1 = runs("H1")
if h1 is not None:
    wHs = [t["p"][3] for t in h1["trace"] if abs(t["p"][0] - h1["zeta"]) < 1e-9
           and t["p"][4] == 0]
    if wHs:
        put("wHmaxLow", max(wHs), "{:.2f}")
        last = max((t for t in h1["trace"] if t["p"][4] == 0 and abs(t["p"][0] - h1["zeta"]) < 1e-9),
                   key=lambda t: t["p"][3])
        put("uLowAtMax", last["u"][0], "{:+.4f}")

# Jacobian check
log = os.path.join(DATA, "check_jacobian.log")
if os.path.exists(log):
    m = re.search(r"worst relative error over all blocks: ([0-9.e+-]+)", open(log).read())
    if m:
        put("jacWorst", float(m.group(1)), "sci0")

# dispersion
dj = os.path.join(DATA, "dispersion.json")
if os.path.exists(dj):
    d = json.load(open(dj))
    k = np.array(d["k"])
    for c in d["cases"]:
        g = np.array(c["maxGrowth"])
        if c["w"] == [1, 1, 1, 0]:
            put("kUnstable", float(k[g > 1e-9].min()), "{:.1f}")
        if c["w"] == [1, 1, 1, 1]:
            put("ampCoeffFull", c["amplitude_k4_coefficient"], "{:.2f}")
        if c["w"] == [1, 1, 1, 0]:
            put("ampCoeffLow", c["amplitude_k4_coefficient"], "{:.2f}")


def tex_name(s):
    # LaTeX macro names: letters only
    s = s.replace("0", "Zero").replace("1", "One").replace("2", "Two").replace("5", "Five")
    return "\\num" + re.sub(r"[^A-Za-z]", "", s)


used = set(re.findall(r"\\num[A-Za-z]+", open(os.path.join(HERE, "report.tex")).read()))
defined = {tex_name(k_) for k_ in macros}
missing = sorted(used - defined)
with open(os.path.join(HERE, "numbers.tex"), "w") as fh:
    fh.write("% generated by collect_numbers.py -- do not edit\n")
    for k_, v in sorted(macros.items()):
        fh.write(f"\\newcommand{{{tex_name(k_)}}}{{{v}}}\n")
    for mname in missing:          # used in report.tex but not (yet) in the data
        fh.write(f"\\newcommand{{{mname}}}{{\\textbf{{??}}}}\n")
if missing:
    print("MISSING (printed as ??):", " ".join(missing))
json.dump(num, open(os.path.join(DATA, "numbers.json"), "w"), indent=1)
for k_, v in sorted(macros.items()):
    print(f"{tex_name(k_):32s} {v}")
