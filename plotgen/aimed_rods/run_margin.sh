#!/usr/bin/env bash
# Regenerates aimed_rods_margin.pdf and aimed_rods_margin.png (the figure of aimed_rods.wls
# with its three panels stacked, for the margin of a book) next to this script.
# The figure is not in the manuscript, so nothing is copied to manuscript/figures/.
# Arguments are passed on (width=<design pt>, default 144; out=<dir>). Runs from any directory.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/aimed_rods_margin.wls" "$@"
