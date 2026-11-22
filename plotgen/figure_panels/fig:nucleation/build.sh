#!/usr/bin/env bash
# Builds nucleation.pdf (fig:nucleation) from its panels and copies it to
# manuscript/figures. Panels: (a) page 2 of plotgen/manual/manual.pdf (drawn by hand),
# (b) plotgen/fluctuating_barier/fluctuating_barrier.pdf (its run.sh). Only assembles;
# the panels are not regenerated. Intermediates go to .tmp/. Runs from any directory.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"
name=nucleation

mkdir -p .tmp
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=.tmp "$name.tex" \
  > /dev/null || { tail -n 30 ".tmp/$name.log"; exit 1; }
grep '^panelrow:' ".tmp/$name.log"
mv ".tmp/$name.pdf" .
cp "$name.pdf" ../../../manuscript/figures/
echo "built $name.pdf, copied to manuscript/figures/$name.pdf"
