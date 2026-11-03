#!/usr/bin/env python3
"""Fast Beamer build: compile only the frames changed since the last build.

Usage (from LaTeX Workshop, cwd = directory of the root file):
    fastbeamer.py [--full] [--context N] <file.tex> <outdir>

--context N also compiles the N frames before and the N frames after each changed frame.

* Not a Beamer file            -> plain latexmk build (same as LaTeX Workshop's default).
* --full, no snapshot yet, preamble changed,
  or a changed input file is used
  in the preamble / by no frame  -> full build; <outdir>/<name>.pdf is the whole deck,
                                  also kept as <outdir>/<name>-full.pdf.
* otherwise                    -> every unchanged frame is blanked out (line numbers
                                  are kept, so SyncTeX still points into <file.tex>),
                                  the rest is compiled, and the result is copied to
                                  <outdir>/<name>.pdf so the viewer shows it.

A frame is compiled when, compared with the last successful build,
  - its text, or the text between it and the previous frame, differs in any byte
    (comments and whitespace included),
  - a neighbouring frame was deleted, or the text after the last frame changed
    (then the frame next to the edit is compiled, so the edit is not lost),
  - an input file it names (figure, \\input file, ...) changed. Input files are the
    project files recorded in latexmk's .fls plus the \\addbibresource files.
"""
import argparse
import difflib
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

LATEXMK = ["latexmk", "-synctex=1", "-interaction=nonstopmode", "-file-line-error", "-pdf"]


def match_brace(text, i):
    """text[i] == '{'; return the index just after the matching '}'."""
    depth = 0
    while i < len(text):
        c = text[i]
        if c == "\\":
            i += 2
            continue
        if c == "%":
            i = text.find("\n", i)
            if i < 0:
                return len(text)
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    return len(text)


TOKEN = re.compile(r"\\begin\{frame\}|\\end\{frame\}|\\frame\s*(?=\{)|(?<!\\)%[^\n]*"
                   r"|\\comment\s*(?=\{)|\\iffalse\b|\\end\{document\}")


def find_frames(text, start):
    """Return [(begin, end, numbered)] spans of top-level frames after `start`."""
    frames, depth, begin = [], 0, None
    pos = start
    while True:
        m = TOKEN.search(text, pos)
        if not m:
            break
        tok = m.group(0)
        if tok.startswith("%"):
            pos = m.end()
        elif tok == r"\end{document}" and depth == 0:  # LaTeX stops here; later frames are dead
            break
        elif tok.startswith(r"\comment"):  # \newcommand{\comment}[1]{}: frames inside are dead
            pos = match_brace(text, m.end())
        elif tok == r"\iffalse":
            f = re.compile(r"\\fi\b").search(text, m.end())
            pos = f.end() if f else len(text)
        elif tok == r"\begin{frame}":
            if depth == 0:
                begin = m.start()
            depth += 1
            pos = m.end()
        elif tok == r"\end{frame}":
            depth -= 1
            pos = m.end()
            if depth == 0 and begin is not None:
                frames.append((begin, pos))
                begin = None
        elif depth == 0:  # \frame{...}
            end = match_brace(text, m.end())
            frames.append((m.start(), end))
            pos = end
        else:
            pos = m.end()
    out = []
    for b, e in frames:
        opts = re.match(r"\\begin\{frame\}\s*(<[^>]*>)?\s*\[([^\]]*)\]", text[b:e])
        out.append((b, e, not (opts and "noframenumbering" in opts.group(2))))
    return out


def split(text):
    """Split a Beamer source into (preamble, [(gap, frame_text)], tail, frame spans)."""
    m = re.search(r"^[^%\n]*?(\\begin\{document\})", text, re.M)
    doc_start = m.start(1)
    frames = find_frames(text, doc_start)
    units, prev = [], doc_start
    for b, e, _ in frames:
        units.append((text[prev:b], text[b:e]))
        prev = e
    return text[:doc_start], units, text[prev:], frames


def run(cmd):
    print("$", " ".join(cmd), flush=True)
    return subprocess.call(cmd)


def latex(tex, outdir, job):
    """Run latexmk; 0 if LaTeX produced the PDF without errors.

    latexmk also fails on BibTeX errors (e.g. a repeated .bib entry) although the PDF is
    fine; those are reported but do not count as a failed build.
    """
    pdf, log = os.path.join(outdir, job + ".pdf"), os.path.join(outdir, job + ".log")
    t0 = os.path.getmtime(pdf) if os.path.exists(pdf) else 0
    rc = run(LATEXMK + [f"-outdir={outdir}", tex])
    if rc == 0:
        return 0
    with open(log, encoding="utf-8", errors="replace") as f:
        latex_error = re.search(r"^!|^\S+\.\w+:\d+: ", f.read(), re.M)
    pdf_new = os.path.exists(pdf) and os.path.getmtime(pdf) > t0
    if latex_error or not pdf_new:
        return rc
    print(f"fastbeamer: latexmk exit code {rc}, but LaTeX ran cleanly (BibTeX error?) -> accepted",
          flush=True)
    return 0


def file_state(path, old=None):
    """(mtime, size, sha1) of a file; sha1 is reused from `old` if mtime and size match."""
    try:
        st = os.stat(path)
    except OSError:
        return None
    if old and old[0] == st.st_mtime_ns and old[1] == st.st_size:
        return old
    with open(path, "rb") as f:
        return [st.st_mtime_ns, st.st_size, hashlib.sha1(f.read()).hexdigest()]


def build_inputs(fls, text, tex, outdir):
    """Project files read by the build (from the .fls) plus the bibliography files."""
    paths = set()
    if os.path.exists(fls):
        with open(fls, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("INPUT "):
                    paths.add(line[6:].strip())
    for m in re.finditer(r"\\(?:addbibresource|bibliography)\{([^}]*)\}", text):
        for b in m.group(1).split(","):
            b = b.strip()
            paths.add(b if b.endswith(".bib") else b + ".bib")
    here, out = os.path.abspath("."), os.path.abspath(outdir)
    rel = set()
    for p in paths:
        a = os.path.abspath(p)
        if a.startswith(here + os.sep) and not a.startswith(out + os.sep) and a != os.path.abspath(tex):
            rel.add(os.path.relpath(a, here))
    return rel


def save_state(tex, snapshot, depfile, inputs, old_deps):
    shutil.copyfile(tex, snapshot)
    deps = {p: file_state(p, old_deps.get(p)) for p in set(old_deps) | inputs}
    with open(depfile, "w") as f:
        json.dump({p: v for p, v in deps.items() if v}, f, indent=0)


def full_build(tex, outdir, name, snapshot, depfile, text, old_deps):
    rc = latex(tex, outdir, name)
    if rc == 0:
        shutil.copyfile(os.path.join(outdir, name + ".pdf"), os.path.join(outdir, name + "-full.pdf"))
        inputs = build_inputs(os.path.join(outdir, name + ".fls"), text, tex, outdir)
        # a full build knows every input: forget files that are no longer used
        save_state(tex, snapshot, depfile, inputs, {p: old_deps[p] for p in inputs & set(old_deps)})
    return rc


def mentions(text, path):
    """Does the source text name this input file (with or without extension)?"""
    stem = os.path.splitext(path)[0]
    keys = {path, stem}
    # animategraphics{..}{media/anim/img}{000}{100} names frames media/anim/img000.png
    base = re.sub(r"\d+$", "", stem)
    if base != stem and not base.endswith("/"):
        keys.add(base)
    return any(k in text for k in keys)


def blank(text):
    """Replace text by its newlines only, so line numbers after it do not move."""
    return "\n" * text.count("\n")


def main():
    ap = argparse.ArgumentParser(description="Compile only the Beamer frames changed since the last build.")
    ap.add_argument("--full", action="store_true", help="compile the whole deck")
    ap.add_argument("--context", type=int, default=0, metavar="N",
                    help="also compile N frames before and after each changed frame")
    ap.add_argument("tex")
    ap.add_argument("outdir")
    a = ap.parse_args()
    force_full, tex, outdir = a.full, a.tex, a.outdir
    name = os.path.splitext(os.path.basename(tex))[0]
    os.makedirs(outdir, exist_ok=True)
    snapshot = os.path.join(outdir, name + ".last.tex")
    depfile = os.path.join(outdir, name + ".deps.json")

    with open(tex, encoding="utf-8") as f:
        text = f.read()

    if not re.search(r"\\documentclass(\[[^\]]*\])?\{beamer\}", text):
        return run(LATEXMK + [f"-outdir={outdir}", tex])

    old_deps = {}
    if os.path.exists(depfile):
        with open(depfile) as f:
            old_deps = json.load(f)

    def full(reason):
        print(f"fastbeamer: {reason} -> full build", flush=True)
        return full_build(tex, outdir, name, snapshot, depfile, text, old_deps)

    if force_full:
        return full("requested")
    if not os.path.exists(snapshot) or not os.path.exists(depfile):
        return full("no previous build")

    with open(snapshot, encoding="utf-8") as f:
        old_text = f.read()
    pre, units, tail, frames = split(text)
    old_pre, old_units, old_tail, _ = split(old_text)
    if pre != old_pre:
        return full("preamble changed")
    if not units:
        return full("no frames")

    keep = set()
    sm = difflib.SequenceMatcher(None, old_units, units, autojunk=False)
    for tag, _, _, j1, j2 in sm.get_opcodes():
        if tag in ("replace", "insert"):
            keep |= set(range(j1, j2))
        elif tag == "delete":  # show the frame that now stands where the deleted one was
            keep.add(min(j1, len(units) - 1))
    if tail != old_tail:
        keep.add(len(units) - 1)

    for path, old in old_deps.items():
        if file_state(path, old) == old:
            continue
        if mentions(pre, path):
            return full(f"{path} changed (used in the preamble)")
        users = [k for k, (gap, body) in enumerate(units) if mentions(gap + body, path)]
        if not users:
            return full(f"{path} changed (not named by any frame)")
        print(f"fastbeamer: {path} changed", flush=True)
        keep |= set(users)

    keep = {j for k in keep for j in range(k - a.context, k + a.context + 1) if 0 <= j < len(units)}
    if not keep:
        print("fastbeamer: nothing changed since the last build", flush=True)
        return 0

    numbers, n = [], 0
    for _, _, numbered in frames:
        n += numbered
        numbers.append(n)
    print("fastbeamer: compiling frame(s)", ", ".join(str(numbers[k]) for k in sorted(keep)), flush=True)

    # Blank unchanged frames; keep everything between frames (sections, groups, settings).
    # The "n / total" footer reads the total from the previous run's .nav; pin it instead.
    doc = re.search(r"^[^%\n]*?(\\begin\{document\})", text, re.M).end(1)
    parts, prev = [text[:doc], r"\makeatletter\gdef\inserttotalframenumber{%d}\makeatother" % n], doc
    for k, (b, e, numbered) in enumerate(frames):
        parts.append(text[prev:b])
        if k in keep:
            # restore the original frame number (setcounter is global, works inside groups)
            parts.append(r"\setcounter{framenumber}{%d}" % (numbers[k] - numbered) + text[b:e])
        else:
            parts.append(blank(text[b:e]))
        prev = e
    # keep the "n / total" footer right: the total is the framenumber at \end{document}
    tail = text[prev:]
    parts.append(re.sub(r"^([^%\n]*?)(\\end\{document\})",
                        lambda m: m.group(1) + r"\setcounter{framenumber}{%d}" % n + m.group(2),
                        tail, count=1, flags=re.M))

    pdir = os.path.join(outdir, "changes")
    os.makedirs(pdir, exist_ok=True)
    ptex = os.path.join(pdir, name + "-changes.tex")
    with open(ptex, "w", encoding="utf-8") as f:
        f.write("".join(parts))
    rc = latex(ptex, pdir, name + "-changes")
    if rc != 0:
        print(f"fastbeamer: build failed, see {os.path.join(pdir, name + '-changes.log')}", flush=True)
        return rc

    shutil.copyfile(os.path.join(pdir, name + "-changes.pdf"), os.path.join(outdir, name + ".pdf"))
    # Point SyncTeX at the real source: same line numbers, different file name.
    src = os.path.join(pdir, name + "-changes.synctex.gz")
    if os.path.exists(src):
        with gzip.open(src, "rt", encoding="utf-8", errors="surrogateescape") as f:
            sync = f.read()
        sync = sync.replace(os.path.abspath(ptex), os.path.abspath(tex))
        sync = re.sub(r"(Input:\d+:)[^\n]*" + re.escape(name + "-changes.tex"),
                      lambda m: m.group(1) + os.path.abspath(tex), sync)
        with gzip.open(os.path.join(outdir, name + ".synctex.gz"), "wt", encoding="utf-8",
                       errors="surrogateescape") as f:
            f.write(sync)
    inputs = build_inputs(os.path.join(pdir, name + "-changes.fls"), text, ptex, outdir)
    save_state(tex, snapshot, depfile, inputs, old_deps)
    return 0


if __name__ == "__main__":
    sys.exit(main())
