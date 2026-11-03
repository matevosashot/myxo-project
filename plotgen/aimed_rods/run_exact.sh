#!/usr/bin/env bash
# Regenerates aimed_rods_exact.pdf and aimed_rods_exact.png (fig:aimed-rods of the manuscript,
# App F: rods aimed at, or side-on to, a point in a uniform nematic and on the head and tail
# sides of a +1/2 defect, on the profile eq:pade-accurate) next to this script, and copies the
# pdf to manuscript/figures/. Arguments are passed on (out=<dir> writes elsewhere and skips
# the copy). Runs from any directory.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/aimed_rods_exact.wls" "$@"

if [[ "$*" != *out=* ]]; then
  cp "$here/aimed_rods_exact.pdf" "$repo/manuscript/figures/aimed_rods_exact.pdf"
  echo "copied aimed_rods_exact.pdf to manuscript/figures/"
fi
