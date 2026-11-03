#!/bin/bash
# LaTeX Workshop recipe: build a latexdiff PDF of <file.tex> against a git commit
# and put it where the viewer looks (<outdir>/<name>.pdf).
#
# Usage (cwd = directory of the root file, as LaTeX Workshop runs it):
#   latexdiff_build.sh <file.tex> <outdir> [commit] [--keep-highlights] [--hide-deleted]
#   commit             commit to compare against, or "staged" for the version in the git
#                      index (staged, not yet committed; equals HEAD when nothing is staged)
#                      (default: $LATEXDIFF_COMMIT, else staged)
#   --keep-highlights  keep \upd, \updb, \updc and updenv in both versions
#                      (default: strip them, so the diff shows only real changes)
#   --hide-deleted     drop the deleted (red) text: only the additions (blue) are shown;
#                      the copy next to the tex file is then diff_<commit>_new.pdf
#
# Build files go to .tmp/diff/ next to the tex file; the diff is also kept as
# diff_<commit>.pdf next to it. Read-only for git: uses only `git show`.
#
# The diff overwrites <outdir>/<name>.pdf, so latexmk's database for <name> is removed:
# the next normal build then recompiles even if the source did not change.
# SyncTeX works: the diff's .synctex.gz is rewritten to point at <file.tex>, with each
# line of the latexdiff output mapped to its line in <file.tex> (line alignment by difflib).

set -euo pipefail

TEX="$1"; OUTDIR="$2"; shift 2
COMMIT="${LATEXDIFF_COMMIT:-staged}"   # set in the recipe's "env" in .vscode/settings.json
STRIP=1
HIDE=0
for a in "$@"; do
  case "$a" in
    --keep-highlights) STRIP=0 ;;
    --hide-deleted) HIDE=1 ;;
    *) COMMIT="$a" ;;
  esac
done

# LaTeX Workshop passes absolute paths: work from the tex file's folder with its bare name
mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)
cd "$(dirname "$TEX")"
TEX=$(basename "$TEX")
NAME=$(basename "$TEX" .tex)
TMP=.tmp/diff
mkdir -p "$TMP" "$OUTDIR"
if [ "$COMMIT" = staged ]; then SHORT=staged; REV=""      # git show :path = index version
else SHORT=$(git rev-parse --short "$COMMIT"); REV="$SHORT"; fi
JOB="diff_${SHORT}"
if [ "$HIDE" = 1 ]; then JOB="${JOB}_new"; fi
REL=$(git rev-parse --show-prefix)   # path of this folder inside the repo

git show "${REV}:${REL}${TEX}" > "$TMP/old.tex"
cp "$TEX" "$TMP/new.tex"

if [ "$STRIP" = 1 ]; then
  # brace-aware removal of the highlight wrappers (content kept)
  python3 - "$TMP/old.tex" "$TMP/new.tex" <<'EOF'
import re, sys
def strip(s):
    for cmd in ('updb', 'updc', 'upd'):
        pat = re.compile(r'\\' + cmd + r'\{')
        out, i = [], 0
        while True:
            m = pat.search(s, i)
            if not m:
                out.append(s[i:]); break
            # skip the definitions in the preamble
            if s[max(0, m.start() - 12):m.start()].endswith(('newcommand{', 'renewcommand{')):
                out.append(s[i:m.end()]); i = m.end(); continue
            out.append(s[i:m.start()])
            j, depth = m.end(), 1
            while depth:
                c = s[j]
                if c == '\\': j += 2; continue
                if c == '{': depth += 1
                elif c == '}': depth -= 1
                j += 1
            out.append(s[m.end():j - 1]); i = j
        s = ''.join(out)
    s = re.sub(r'\\begin\{updenv\}\s*\n?', '', s)
    s = re.sub(r'\\end\{updenv\}\s*\n?', '', s)
    return s
for f in sys.argv[1:]:
    t = open(f).read()
    while True:          # repeat: tags can be nested
        u = strip(t)
        if u == t: break
        t = u
    open(f, 'w').write(t)
EOF
fi

# markup: colors only (added blue, deleted red), no underline or strike-through;
# deleted text is also smaller (\footnotesize; not in math, where size commands fail)
cat > "$TMP/diffpreamble.tex" <<'EOF'
\RequirePackage{color}
\definecolor{diffadd}{rgb}{0,0.2,0.85}
\definecolor{diffdel}{rgb}{0.85,0,0}
\providecommand{\DIFadd}[1]{{\protect\color{diffadd}#1}}
\DeclareRobustCommand{\DIFdelsize}{\ifmmode\else\footnotesize\fi}
\providecommand{\DIFdel}[1]{{\protect\color{diffdel}\DIFdelsize#1}}
\providecommand{\DIFaddbegin}{}
\providecommand{\DIFaddend}{}
\providecommand{\DIFdelbegin}{}
\providecommand{\DIFdelend}{}
\providecommand{\DIFmodbegin}{}
\providecommand{\DIFmodend}{}
\providecommand{\DIFaddFL}[1]{\DIFadd{#1}}
\providecommand{\DIFdelFL}[1]{\DIFdel{#1}}
\providecommand{\DIFaddbeginFL}{}
\providecommand{\DIFaddendFL}{}
\providecommand{\DIFdelbeginFL}{}
\providecommand{\DIFdelendFL}{}
EOF
if [ "$HIDE" = 1 ]; then   # deleted text is swallowed
  sed -i 's/^\\providecommand{\\DIFdel}\[1\].*/\\providecommand{\\DIFdel}[1]{}/' "$TMP/diffpreamble.tex"
fi
echo "latexdiff_build: $TEX vs $SHORT"
latexdiff --math-markup=coarse --preamble="$TMP/diffpreamble.tex" \
  "$TMP/old.tex" "$TMP/new.tex" > "$TMP/$JOB.tex" 2> "$TMP/latexdiff.err"
if [ "$HIDE" = 1 ]; then
  # cut whole deleted blocks (also deleted equations and figures, which would leave
  # empty displays); their newlines stay, commented so that no blank line (= \par)
  # appears, and lines still align for SyncTeX
  python3 - "$TMP/$JOB.tex" <<'EOF'
import re, sys
f = sys.argv[1]
t = open(f).read()
pre, sep, body = t.partition('\\begin{document}')
body = re.sub(r'\\DIFdelbegin(FL)?(?![A-Za-z]).*?\\DIFdelend(FL)?(?![A-Za-z])',
              lambda m: '%\n' * m.group(0).count('\n'), body, flags=re.S)
open(f, 'w').write(pre + sep + body)
EOF
fi

# compile from this folder so figures and the .bib resolve; drop stale aux files first.
# -f: go on through LaTeX errors, otherwise one error stops latexmk after the first pass
# and, with the aux files fresh, every reference and citation is left as "?"
rm -f "$TMP/$JOB".{aux,toc,out,fls,fdb_latexmk,bbl,blg,log}
latexmk -pdf -f -synctex=1 -outdir="$TMP" -interaction=nonstopmode -file-line-error "$TMP/$JOB.tex" || true
if [ ! -f "$TMP/$JOB.pdf" ]; then
  echo "latexdiff_build: build failed; see $TMP/$JOB.log" >&2; exit 1
fi
echo "latexdiff_build: LaTeX errors: $(grep -c ':[0-9]*: \|^! ' "$TMP/$JOB.log" || true)"

cp "$TMP/$JOB.pdf" "./$JOB.pdf"
cp "$TMP/$JOB.pdf" "$OUTDIR/$NAME.pdf"
# the viewer's PDF is now the diff: make the next normal build recompile
rm -f "$OUTDIR/$NAME.fdb_latexmk"

# SyncTeX for the diff, pointing at the real source
python3 - "$TMP/$JOB.synctex.gz" "$TMP/$JOB.tex" "$(pwd)/$TEX" "$OUTDIR/$NAME.synctex.gz" <<'EOF'
import difflib, gzip, os, re, sys
syn, difftex, src, dest = sys.argv[1:]
d = open(difftex, errors='replace').read().splitlines()
s = open(src, errors='replace').read().splitlines()
# line map diff -> source (1-based); unmatched lines go proportionally through
# their replaced block, or to the source line where an insertion sits
lmap = {}
sm = difflib.SequenceMatcher(None, d, s, autojunk=False)
for op, i1, i2, j1, j2 in sm.get_opcodes():
    for k in range(i1, i2):
        if op == 'equal':
            j = j1 + (k - i1)
        elif op == 'replace':
            j = j1 + (k - i1) * (j2 - j1) // (i2 - i1)
        else:   # 'delete': diff-only line (old text), no source counterpart
            j = min(j1, len(s) - 1)
        lmap[k + 1] = j + 1
with gzip.open(syn, 'rt', errors='replace') as f:
    lines = f.read().split('\n')
tags = set()
for n, l in enumerate(lines):
    m = re.match(r'Input:(\d+):(.*)$', l)
    if m and os.path.basename(m.group(2)) == os.path.basename(difftex):
        tags.add(m.group(1))
        lines[n] = f'Input:{m.group(1)}:{src}'
rec = re.compile(r'^([\[(vhxkg$])(\d+),(\d+)')
def fix(m):
    if m.group(2) not in tags:
        return m.group(0)
    return f'{m.group(1)}{m.group(2)},{lmap.get(int(m.group(3)), int(m.group(3)))}'
lines = [rec.sub(fix, l) for l in lines]
with gzip.open(dest, 'wt') as f:
    f.write('\n'.join(lines))
EOF
echo "latexdiff_build: written $(pwd)/$JOB.pdf and $OUTDIR/$NAME.pdf"
