#!/usr/bin/env python3
"""Parse the batch log produced by run_many.sh into tidy CSV tables.

Input : ../out_many.txt  (or a path given on the command line)
Output: runs.csv           one row per attempted run (incl. failures)
        timings_long.csv   one row per (run, stage) timing-table entry
        timings_wide.csv   one row per run, one column per stage

The log is a concatenation of run blocks.  Each block after the first is
introduced by

    === fbox=20, Nint=40 -- Wed Aug 12 06:09:21 PM CEST 2026 ===

and ends either with a timing table

    === timing breakdown (Nint = 40, n = 4800) ===
    stage                                                s   %tracked
    ------------------------------------------------------------------
    [6] Eigensystem (dense, nonsymmetric)           123.45    74.8  <== >4 s
    ...
    ------------------------------------------------------------------
    TOTAL tracked                                   165.00   100.0

or with

    === FAILED (exit 131): fbox=20, Nint=120 ===

The first block in the file has no header line (it predates the batch loop);
its parameters are recovered from the ``SCRIPT mode;`` line that every run
prints.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_LOG = HERE.parent / "out_many.txt"

RE_HEADER = re.compile(r"^=== fbox=(\d+), Nint=(\d+) -- (.*?) ===\s*$")
RE_FAILED = re.compile(r"^=== FAILED \(exit (\d+)\): fbox=(\d+), Nint=(\d+) ===\s*$")
RE_SCRIPT = re.compile(r"SCRIPT mode; fbox=(\d+), Nint=(\d+), plotBox=(\d+), saveData=(\w+)")
RE_BREAKDOWN = re.compile(r"^=== timing breakdown \(Nint = (\d+), n = (\d+)\) ===\s*$")
RE_TOTAL = re.compile(r"^TOTAL tracked\s+([\d.]+)\s+([\d.]+)\s*$")
RE_GRID = re.compile(r"^\s*\[2\] Grid: fbox=(\d+), Nint=(\d+), h=([\d.]+), state-dim=(\d+)")
RE_SAVED = re.compile(
    r"^Saved: (\S+?)(?:\s+\(([\d.]+) MiB, ([\d.]+) s\))?\s*$"
)
RE_STAGE = re.compile(r"^(\[\d\] .+?)\s\s+([\d.]+)\s+([\d.]+)\s*(<==.*)?$")
RE_SOLVE_TOTAL = re.compile(r"^\s*total solve time = ([\d.eE*^+-]+) s")
RE_MAXRE = re.compile(r"^\s*max Re\(.\) = (\S+)")


def mma_float(text: str) -> float:
    """Convert a Mathematica-formatted real (1.23*^-9) to a python float."""
    return float(text.replace("*^", "e"))


def split_blocks(lines):
    """Yield (header_info, block_lines) for every run block in the log.

    header_info is a dict with fbox/Nint/timestamp, or None for the leading
    headerless block.
    """
    current_header = None
    current: list[str] = []
    for line in lines:
        m = RE_HEADER.match(line)
        if m:
            if current:
                yield current_header, current
            current_header = {
                "fbox": int(m.group(1)),
                "Nint": int(m.group(2)),
                "timestamp": m.group(3),
            }
            current = []
        else:
            current.append(line)
    if current:
        yield current_header, current


def parse_block(header, lines):
    """Turn one run block into (run_record, [stage_records])."""
    run = {
        "fbox": header["fbox"] if header else None,
        "Nint": header["Nint"] if header else None,
        "timestamp": header["timestamp"] if header else "",
        "n": None,
        "h": None,
        "status": "ok",
        "exit_code": "",
        "total_tracked_s": None,
        "lyapunov_solve_s": None,
        "max_re_lambda": None,
        "cov_size_MiB": None,
        "cov_save_s": None,
        "cov_path": "",
    }
    stages = []
    in_table = False

    for line in lines:
        # ---- run parameters -------------------------------------------------
        m = RE_SCRIPT.search(line)
        if m:
            run["fbox"] = int(m.group(1))
            run["Nint"] = int(m.group(2))
            continue

        m = RE_GRID.match(line)
        if m:
            run["h"] = float(m.group(3))
            run["n"] = int(m.group(4))
            continue

        m = RE_SOLVE_TOTAL.match(line)
        if m:
            run["lyapunov_solve_s"] = mma_float(m.group(1))
            continue

        m = RE_MAXRE.match(line)
        if m:
            run["max_re_lambda"] = mma_float(m.group(1))
            continue

        # ---- saved artefacts ------------------------------------------------
        m = RE_SAVED.match(line)
        if m and m.group(1).endswith(".bin"):
            run["cov_path"] = m.group(1)
            if m.group(2):
                run["cov_size_MiB"] = float(m.group(2))
                run["cov_save_s"] = float(m.group(3))
            continue

        # ---- failure marker -------------------------------------------------
        m = RE_FAILED.match(line)
        if m:
            run["status"] = "failed"
            run["exit_code"] = int(m.group(1))
            continue

        # ---- timing table ---------------------------------------------------
        m = RE_BREAKDOWN.match(line)
        if m:
            in_table = True
            run["n"] = int(m.group(2))
            continue

        if not in_table:
            continue

        m = RE_TOTAL.match(line)
        if m:
            run["total_tracked_s"] = float(m.group(1))
            in_table = False
            continue

        m = RE_STAGE.match(line)
        if m:
            stage = m.group(1).strip()
            stages.append(
                {
                    "stage_group": stage[1],           # the k in "[k] ..."
                    "stage": stage,
                    "seconds": float(m.group(2)),
                    "pct_tracked": float(m.group(3)),
                }
            )

    if run["status"] == "ok" and run["total_tracked_s"] is None:
        run["status"] = "incomplete"
    return run, stages


def main(argv):
    log_path = Path(argv[1]) if len(argv) > 1 else DEFAULT_LOG
    lines = log_path.read_text(errors="replace").splitlines()

    runs, long_rows = [], []
    for header, block in split_blocks(lines):
        run, stages = parse_block(header, block)
        if run["Nint"] is None and run["status"] == "ok":
            continue                                    # not a run block at all
        run_id = f"box{run['fbox']}_N{run['Nint']}"
        run["run_id"] = run_id
        runs.append(run)
        for s in stages:
            long_rows.append(
                {
                    "run_id": run_id,
                    "fbox": run["fbox"],
                    "Nint": run["Nint"],
                    "n": run["n"],
                    **s,
                }
            )

    runs.sort(key=lambda r: (r["fbox"], r["Nint"]))

    # ---- runs.csv -----------------------------------------------------------
    run_cols = [
        "run_id", "fbox", "Nint", "n", "h", "status", "exit_code",
        "total_tracked_s", "lyapunov_solve_s", "max_re_lambda",
        "cov_size_MiB", "cov_save_s", "cov_path", "timestamp",
    ]
    write_csv(HERE / "runs.csv", run_cols, runs)

    # ---- timings_long.csv ---------------------------------------------------
    long_cols = ["run_id", "fbox", "Nint", "n", "stage_group", "stage",
                 "seconds", "pct_tracked"]
    long_rows.sort(key=lambda r: (r["fbox"], r["Nint"], -r["seconds"]))
    write_csv(HERE / "timings_long.csv", long_cols, long_rows)

    # ---- timings_wide.csv ---------------------------------------------------
    # Stage columns in a stable order: by the [k] tag, then by the largest time
    # the stage ever reaches (so the expensive stages come first within a tag).
    peak = {}
    for r in long_rows:
        peak[r["stage"]] = max(peak.get(r["stage"], 0.0), r["seconds"])
    stage_names = sorted(peak, key=lambda s: (s[1], -peak[s]))

    wide_rows = []
    for run in runs:
        row = {k: run[k] for k in
               ("run_id", "fbox", "Nint", "n", "status",
                "total_tracked_s", "cov_size_MiB")}
        times = {r["stage"]: r["seconds"]
                 for r in long_rows if r["run_id"] == run["run_id"]}
        for s in stage_names:
            row[s] = times.get(s, "")
        wide_rows.append(row)
    wide_cols = ["run_id", "fbox", "Nint", "n", "status",
                 "total_tracked_s", "cov_size_MiB"] + stage_names
    write_csv(HERE / "timings_wide.csv", wide_cols, wide_rows)

    ok = [r for r in runs if r["status"] == "ok"]
    print(f"parsed {len(runs)} runs ({len(ok)} complete, "
          f"{len(runs) - len(ok)} failed/incomplete), "
          f"{len(long_rows)} stage rows, {len(stage_names)} distinct stages")
    for r in runs:
        if r["status"] != "ok":
            print(f"  ! {r['run_id']}: {r['status']} (exit {r['exit_code']})")


def write_csv(path, columns, rows):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path.name}  ({len(rows)} rows)")


if __name__ == "__main__":
    main(sys.argv)
