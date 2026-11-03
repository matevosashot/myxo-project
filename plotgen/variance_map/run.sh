#!/usr/bin/env bash
# Regenerates variance_map.pdf and variance_map.png (fig:variance3d: the density variance
# around a +1/2 defect as a surface and as two maps with the circle of eq:p-layer-radial, on
# one colour scale) next to this script, and copies the pdf to manuscript/figures/.
# Arguments are passed on: draft, analytic, r=<radius in um>, angle=<deg>, out=<dir>.
# A draft writes variance_map_draft.png only and copies nothing. Runs from any directory.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen
export GALLIUM_DRIVER=softpipe      # llvmpipe crashes the front end on 3D rendering

wolframscript -file "$here/variance_map.wls" "$@"

draft=0; custom_out=0
for a in "$@"; do
  [[ $a == draft ]] && draft=1
  [[ $a == out=* ]] && custom_out=1
done
if [[ $draft == 0 && $custom_out == 0 ]]; then
  cp "$here/variance_map.pdf" "$here/../../manuscript/figures/variance_map.pdf"
  echo "copied variance_map.pdf to manuscript/figures/"
fi
