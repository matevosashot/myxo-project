#!/usr/bin/env bash
# Builds staggered_grid.pdf (fig:staggered-grid) next to this script and copies it to
# manuscript/figures/. Runs from any directory.
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"
pdflatex -interaction=nonstopmode -halt-on-error staggered_grid.tex >/dev/null
rm -f staggered_grid.aux staggered_grid.log
cp staggered_grid.pdf "$here/../../manuscript/figures/staggered_grid.pdf"
echo "wrote $here/staggered_grid.pdf and manuscript/figures/staggered_grid.pdf"
