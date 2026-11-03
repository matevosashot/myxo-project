#!/usr/bin/env bash
# Regenerates layering_profile.pdf and layering_profile.png (fig:defects-layering(c): the
# circle-averaged layering probability about r_max and about the core, with the points of
# Copenhagen et al. 2021 Fig. 3f) next to this script, from the N_int = 200 ssDir Lyapunov run.
# Without arguments it also copies the pdf to manuscript/figures/final.pdf, which the author
# places by hand into figures/handmade/main.pdf (page 2).
# Arguments are passed on: N=<n>, pcrit=<Pa um>, norm=pedestal, out=<dir> (see
# layering_profile.wls); with any argument nothing is copied.
# Runs from any directory; about 30 s.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/layering_profile.wls" "$@"

if [[ $# == 0 ]]; then
  cp "$here/layering_profile.pdf" "$here/../../manuscript/figures/final.pdf"
  echo "copied layering_profile.pdf to manuscript/figures/final.pdf"
fi
