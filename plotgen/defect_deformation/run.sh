#!/usr/bin/env bash
# Regenerates strength, defectangle and activeforce (.pdf, .png; the three panels of
# fig:comoving-deformation) next to this script, or in out=<dir>. Copies nothing into
# manuscript/figures/. Arguments are passed on: repro, out=<dir> (see defect_deformation.wls).
# Runs from any directory; about a minute.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/defect_deformation.wls" "$@"
