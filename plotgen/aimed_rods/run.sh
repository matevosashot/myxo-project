#!/usr/bin/env bash
# Regenerates aimed_rods.pdf and aimed_rods.png (the physical picture of the head/tail
# asymmetry of the density variance: rods aimed at, or side-on to, a spot in a uniform
# nematic and on the head and tail sides of a +1/2 defect) next to this script.
# The figure is not in the manuscript, so nothing is copied to manuscript/figures/.
# Arguments are passed on (out=<dir> writes elsewhere). Runs from any directory.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/aimed_rods.wls" "$@"
