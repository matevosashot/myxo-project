#!/usr/bin/env bash
# Regenerates circle_average.pdf and circle_average.png (the circle of eq:p-layer-radial on
# density plots of the variance of fig:variance3d) next to this script.
# Arguments are passed on: draft, analytic, r=<radius in um>, angle=<deg>, out=<dir>.
# Runs from any directory.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/circle_average.wls" "$@"
