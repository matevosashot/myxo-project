#!/usr/bin/env bash
# Regenerates varrho.pdf and varrho.png (fig:variance: the density variance over its pedestal
# along and across the defect axis -- analytic at delta Q = 0 on the static profile, the same
# formula on the numerical steady state of the moving defect, and the Lyapunov solution) next
# to this script. Copies nothing into manuscript/figures/.
# Arguments are passed on: check, ssDir, nonew, Rc=<um>, out=<dir> (see varrho.wls).
# Runs from any directory; about 1 min (about 3 min with check).
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/varrho.wls" "$@"
