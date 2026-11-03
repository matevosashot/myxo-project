#!/usr/bin/env bash
# Regenerates drift_deformation.pdf and drift_deformation.png (fig:drift-deformation of
# app:drift: director rotation near the core without and with the closure term T) next to
# this script, and copies the pdf to manuscript/figures/. Arguments are passed on
# (out=<dir> writes elsewhere and skips the copy). Runs from any directory.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/drift_deformation.wls" "$@"
if [ $# -eq 0 ]; then
  cp "$here/drift_deformation.pdf" "$here/../../manuscript/figures/"
  echo "copied drift_deformation.pdf to manuscript/figures/"
fi
