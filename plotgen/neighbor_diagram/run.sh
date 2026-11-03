#!/usr/bin/env bash
# Builds the three alternative neighbor diagrams for eq:total-gauss / eq:omega-weight
# (neighbor_diagram_v1/v2/v3.pdf and 300 dpi PNG renders) next to this script.
# omega_data.py first writes the data of v2 (gray rings) and v3 (the curve 2 pi d Omega).
# Intermediate files go to .tmp/. v3 is fig:neighbors of the manuscript and is copied to
# manuscript/figures/neighbor_diagram_v3.pdf. Runs from any directory; pass version numbers
# to build only some, e.g. ./run.sh 1 3.
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"
mkdir -p .tmp

python3 omega_data.py

versions=("$@")
[[ ${#versions[@]} -eq 0 ]] && versions=(1 2 3)
for v in "${versions[@]}"; do
  name="neighbor_diagram_v$v"
  if ! pdflatex -interaction=nonstopmode -halt-on-error -output-directory=.tmp "$name.tex" \
      >".tmp/$name.stdout"; then
    grep -A5 '^!' ".tmp/$name.log" || tail -20 ".tmp/$name.stdout"
    exit 1
  fi
  mv ".tmp/$name.pdf" "$name.pdf"
  # poppler pdftoppm (the xpdf one in /usr/local/bin has no -png)
  /usr/bin/pdftoppm -r 300 -png -singlefile "$name.pdf" "$name"
  echo "wrote $here/$name.pdf and $name.png"
  if [[ $v == 3 ]]; then
    cp "$name.pdf" ../../manuscript/figures/
    echo "copied $name.pdf to manuscript/figures/"
  fi
done
