#!/usr/bin/env python3
"""Log-log timing scaling of the Lyapunov sweep, from the runs' own metadata.

    python3 lyapunov_solver/plot_scaling.py                       # replot from cache
    python3 lyapunov_solver/plot_scaling.py --refresh              # re-read the meta files
    python3 lyapunov_solver/plot_scaling.py --refresh --jobs 50:1208115,100:1208116

WHERE THE NUMBERS COME FROM -- read this before trusting a point
----------------------------------------------------------------
The solve time is read from `covariance_<tag>_meta.m`, NOT from the SLURM logs.
That matters: for a long solve the task log goes silent.  Every Nint=180 log in
the first sweep stops dead after "[w] write A,Q" and never prints the timing
breakdown -- the kernel sits in one blocking RunProcess for ~12.8 h and its
buffered stdout never reaches the file, though it keeps computing and writes
every output correctly.  The meta file is written by the solver itself and is
unaffected, so it is the reliable source at every size.  Expect the same silence
at Nint >= 200 (~24 h and ~46 h of no output).

The end-to-end wall time comes from `sacct` and needs the array job ids, so it is
optional: without --jobs (or a cached copy) only the solve series is drawn.

WHAT IT SHOWS
-------------
Measured on the 2026-09 sweep: the solve scales as Nint^5.4 = n^2.69, short of
the ideal Nint^6 = n^3.  The exponent is not wrong -- BLAS efficiency rises with
size, so the constant falls across this range; the asymptote is still n^3.  The
solve is 48% of the task at Nint=50 but 93-98% from Nint=100 up, so above ~100
the task time IS the solve time.
"""

import argparse
import collections
import json
import os
import statistics as st
import subprocess
import sys

DEFAULT_DATA = ("/data/biophys/ashmat/data/myxo-project/"
                "fluctuations-set2_2_comove_optimized")
HERE = os.path.dirname(os.path.abspath(__file__))
WOLFRAM = "/usr/local/math/math1410/Executables/wolframscript"
SOLVE_KEY = "[6] python schur + recursive Bartels-Stewart"

# dataviz default categorical slots 1-2, validated for CVD separation
# (normal dE 33.6, protan 26.5, deutan 33.3, tritan 28.7; floors 15 / 8).
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#3d3d3d", "#6b6b6b", "#e4e4e1", "#fcfcfb"


def extract_solve(data_dir, cache):
    """Pull the solve timing out of every covariance_*_meta.m via wolframscript."""
    # NB: never name a Mathematica variable `D` or `N` here -- both are Protected
    # builtins and the assignment fails with a message that points somewhere else.
    code = r'''
    dir = "%s";
    fs = FileNames["covariance_*_meta.m", dir];
    rows = Table[Module[{d, t, s},
       d = Quiet@Check[Import[f], $Failed];
       If[d === $Failed || ! KeyExistsQ[d, "timings"], Nothing,
        t = Association[Rule @@@ d["timings"]];
        s = Lookup[t, "%s", Missing[]];
        If[NumericQ[s],
         <|"Nint" -> d["Nint"], "n" -> Lookup[d, "n", Missing[]], "solve" -> s,
           "tag" -> Lookup[d, "runTag", ""]|>, Nothing]]], {f, fs}];
    Export["%s", rows, "JSON"];
    Print["extracted ", Length[rows], " of ", Length[fs], " meta files"];
    ''' % (data_dir, SOLVE_KEY, cache)
    wolf = WOLFRAM if os.path.exists(WOLFRAM) else "wolframscript"
    r = subprocess.run([wolf, "-code", code], capture_output=True, text=True)
    # `wolframscript -code` echoes the value of the last expression, and Print
    # returns Null, so a bare "Null" line always trails the real output.
    sys.stdout.write("\n".join(l for l in r.stdout.splitlines()
                                if l.strip() != "Null") + "\n")
    if not os.path.exists(cache):
        sys.exit("wolframscript produced no cache; stderr:\n" + r.stderr[:2000])
    return json.load(open(cache))


def extract_wall(jobs, cache):
    """Median end-to-end task wall from sacct, per Nint. jobs = {Nint: jobid}."""
    def secs(e):
        d = 0
        if "-" in e:
            d, e = e.split("-")
            d = int(d)
        h, m, s = e.split(":")
        return d * 86400 + int(h) * 3600 + int(m) * 60 + float(s)

    out = collections.defaultdict(list)
    for N, jid in jobs.items():
        r = subprocess.run(["sacct", "-j", str(jid), "-n", "-P",
                            "-o", "JobID,State,Elapsed"],
                           capture_output=True, text=True).stdout
        for line in r.splitlines():
            p = line.split("|")
            if len(p) < 3 or ".batch" in p[0] or ".extern" in p[0]:
                continue
            if p[1] != "COMPLETED":
                continue
            out[int(N)].append(secs(p[2]))
    json.dump({str(k): v for k, v in out.items()}, open(cache, "w"))
    return {int(k): v for k, v in out.items()}


def slope(d, lo=None):
    import numpy as np
    Ns = sorted(k for k in d if lo is None or k >= lo)
    if len(Ns) < 2:
        return float("nan")
    return float(np.polyfit(np.log([float(n) for n in Ns]),
                            np.log([st.median(d[n]) for n in Ns]), 1)[0])


def plot(solve, wall, out_path):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p_s, p_s_hi = slope(solve), slope(solve, 100)
    series = [(wall, ORANGE, "End-to-end task wall")] if wall else []
    series.append((solve, BLUE, "Lyapunov solve only"))

    fig = plt.figure(figsize=(9.0, 6.9), dpi=170)
    fig.patch.set_facecolor(SURF)
    ax = fig.add_axes([0.105, 0.115, 0.72, 0.74])
    ax.set_facecolor(SURF)

    # individual tasks behind the medians -- the spread is real (node variation)
    for d, c, _ in series:
        for n, v in d.items():
            ax.scatter([n] * len(v), v, s=12, color=c, alpha=0.20,
                       linewidths=0, zorder=2)
    for d, c, lab in series:
        Ns = sorted(d)
        ax.plot(Ns, [st.median(d[n]) for n in Ns], "-o", color=c, lw=2, ms=8,
                mfc=c, mec=SURF, mew=2, label=lab, zorder=4)

    ref = wall or solve
    anc_x = max(ref)
    anc_y = st.median(ref[anc_x])
    lo, hi = min(min(solve), min(ref)), anc_x
    xa = np.array([lo * 0.96, hi * 1.35])
    ax.plot(xa, anc_y * (xa / anc_x) ** 6, "--", color=MUTED, lw=1.4, zorder=1)
    mid = (lo * hi) ** 0.5
    ax.annotate(r"ideal  $O(n^3)=N_{\rm int}^{6}$",
                xy=(mid, anc_y * (mid / anc_x) ** 6), xytext=(6, -30),
                textcoords="offset points", color=MUTED, fontsize=10.5,
                ha="left", va="top")

    # direct labels at the right end, in ink not series colour
    for d, _, lab in series:
        n = max(d)
        dy = 10 if lab.startswith("End") else -16
        ax.annotate(lab.replace(" task", "\ntask").replace(" solve", "\nsolve"),
                    xy=(n, st.median(d[n])), xytext=(11, dy),
                    textcoords="offset points", color=INK2, fontsize=10.5,
                    ha="left", va="center")

    ax.set_xscale("log")
    ax.set_yscale("log")
    Ns_all = sorted(set(solve) | set(wall or {}))
    major = [n for n in Ns_all if n in (50, 100, 150, 180, 200, 230)] or Ns_all
    minor = [n for n in Ns_all if n not in major]
    ax.set_xticks(major)
    ax.set_xticklabels([str(n) for n in major])
    ax.set_xticks(minor, minor=True)
    ax.set_xticklabels([], minor=True)
    ax.set_yticks([60, 300, 1800, 3600, 3 * 3600, 12 * 3600, 48 * 3600])
    ax.set_yticklabels(["1 min", "5 min", "30 min", "1 h", "3 h", "12 h", "2 d"])
    ax.set_yticks([], minor=True)
    ax.set_xlim(min(Ns_all) * 0.92, max(Ns_all) * 1.29)
    ax.set_xlabel(r"$N_{\rm int}$    (state dimension $n=3N_{\rm int}^2$)",
                  color=INK2, fontsize=11.5, labelpad=8)
    ax.set_ylabel("median time per task", color=INK2, fontsize=11.5, labelpad=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=10.5, length=0)
    ax.grid(True, which="major", color=GRID, lw=0.9)
    ax.grid(True, which="minor", axis="x", color=GRID, lw=0.6, alpha=0.6)
    ax.set_axisbelow(True)
    if len(series) > 1:
        ax.legend(loc="upper left", frameon=False, fontsize=10.5,
                  labelcolor=INK2, handlelength=1.8, borderpad=0.2)

    ntask = sum(len(v) for v in solve.values())
    fig.text(0.105, 0.955, "Lyapunov sweep: time scaling with grid size",
             color=INK, fontsize=15, weight="bold", ha="left", va="top")
    fig.text(0.105, 0.905,
             f"solve scales as $t\\propto N_{{\\rm int}}^{{{p_s:.1f}}}$ "
             f"($n^{{{p_s / 2:.2f}}}$), $N_{{\\rm int}}^{{{p_s_hi:.1f}}}$ for "
             f"$N_{{\\rm int}}\\geq100$ — short of the ideal $N_{{\\rm int}}^{{6}}$"
             f"\n{ntask} tasks; faint dots are individual runs"
             + (f" (unlabelled minor ticks: {', '.join(map(str, minor))})" if minor else ""),
             color=MUTED, fontsize=10.5, ha="left", va="top")
    fig.savefig(out_path, facecolor=SURF)
    return p_s, p_s_hi


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=DEFAULT_DATA,
                    help="directory holding covariance_*_meta.m")
    ap.add_argument("--jobs", default="",
                    help="Nint:jobid pairs for sacct wall times, e.g. 50:1208115,100:1208116")
    ap.add_argument("--out", default=os.path.join(HERE, "nint_scaling.png"))
    ap.add_argument("--cache-dir", default=HERE)
    ap.add_argument("--refresh", action="store_true",
                    help="re-read the meta files (and sacct) instead of using the cache")
    a = ap.parse_args()

    solve_cache = os.path.join(a.cache_dir, "scaling_solve.json")
    wall_cache = os.path.join(a.cache_dir, "scaling_wall.json")

    if a.refresh or not os.path.exists(solve_cache):
        rows = extract_solve(a.data_dir, solve_cache)
    else:
        rows = json.load(open(solve_cache))
    solve = collections.defaultdict(list)
    for r in rows:
        solve[int(r["Nint"])].append(float(r["solve"]))
    if not solve:
        sys.exit("no solve timings found — check --data-dir")

    wall = {}
    if a.jobs and (a.refresh or not os.path.exists(wall_cache)):
        jobs = dict(p.split(":") for p in a.jobs.split(","))
        wall = extract_wall(jobs, wall_cache)
    elif os.path.exists(wall_cache):
        wall = {int(k): v for k, v in json.load(open(wall_cache)).items()}

    p_s, p_s_hi = plot(solve, wall, a.out)

    print(f"\nsolve: Nint^{p_s:.2f} (n^{p_s / 2:.2f})   "
          f"Nint>=100: Nint^{p_s_hi:.2f} (n^{p_s_hi / 2:.2f})\n")
    print(f"{'Nint':>5} {'n':>7} {'runs':>5} {'solve s':>10} {'h':>7} {'wall s':>10} {'solve/wall':>11}")
    for n in sorted(solve):
        ms = st.median(solve[n])
        w = st.median(wall[n]) if n in wall else None
        print(f"{n:>5} {3 * n * n:>7} {len(solve[n]):>5} {ms:>10.0f} {ms / 3600:>7.2f} "
              f"{(f'{w:.0f}' if w else '-'):>10} {(f'{ms / w:.0%}' if w else '-'):>11}")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
