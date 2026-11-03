#!/usr/bin/env bash
# Regenerates omega_weight.pdf and omega_weight.png (the weight Omega of eq:variance-neighbors,
# app:variance-neighbors: the radial weight 2 pi w Omega(w) against w, single column) next to
# this script, and copies the pdf to manuscript/figures/. Arguments are passed on (out=<dir>
# writes elsewhere and skips the copy). Runs from any directory.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/omega_weight.wls" "$@"

if [[ "$*" != *out=* ]]; then
  cp "$here/omega_weight.pdf" "$repo/manuscript/figures/omega_weight.pdf"
  echo "copied omega_weight.pdf to manuscript/figures/"
fi
