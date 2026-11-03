#!/usr/bin/env bash
# Regenerates rho1d_Dirichlet, rho1d_Dirichlet_Neumann and nematicaxis (.pdf, .png, _high.png;
# fig:rho1d-dirichlet, fig:rho1d-dirichlet-neumann, fig:nematicaxis) next to this script, or in
# out=<dir>. Copies nothing into manuscript/figures/. Arguments are passed on: resolve,
# out=<dir> (see rho1d.wls). Runs from any directory; a few minutes for the solves, seconds
# with the cache.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/rho1d.wls" "$@"
