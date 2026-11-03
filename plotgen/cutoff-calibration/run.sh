#!/usr/bin/env bash
# Regenerates cutoff_calibration.{pdf,png} (fig:cutoff-calibration, l_eff/l_n against dx/l_n)
# next to this script and copies them to manuscript/figures/; out=<dir> writes there instead
# and copies nothing. Runs from any directory.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/cutoff_calibration.wls" "$@"
if [ $# -eq 0 ]; then              # out=<dir> writes elsewhere and skips the copy
  cp "$here/cutoff_calibration.pdf" "$here/cutoff_calibration.png" "$here/../../manuscript/figures/"
  echo "copied cutoff_calibration.{pdf,png} to manuscript/figures/"
fi
