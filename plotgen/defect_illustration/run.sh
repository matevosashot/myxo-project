#!/usr/bin/env bash
# Regenerates the +1/2 and -1/2 defect icons (defect_plus/minus.{pdf,png}) next to this
# script. Runs from any directory.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen
export GALLIUM_DRIVER=softpipe      # llvmpipe crashes the front end in Rasterize

wolframscript -file "$here/defect_illustration.wls"
