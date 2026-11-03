#!/usr/bin/env bash
# Regenerates the figures of the fluctuating active force (app:lyapunov-results) next to this
# script: fActiveMeanAndStd.pdf (fig:force-axis) and f2.png / f2_high.png (fig:force-map), plus
# fxfx/fyfy pngs, from the ssDir N_int = 200 Lyapunov run; then copies fActiveMeanAndStd.pdf,
# f2.png and f2_high.png to manuscript/figures/.
# Arguments are passed on: N=<Nint>, meanSS, cache=<file.mx>, out=<dir> (see force_fluct.wls).
# With out= or N= nothing is copied. Runs from any directory; about 8 min from a cold read of
# the 115 GB covariance (0.4 GB of memory), about 1 min with cache=.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen
export GALLIUM_DRIVER=softpipe      # llvmpipe crashes the front end on raster export

wolframscript -file "$here/force_fluct.wls" "$@"

copy=1
for a in "$@"; do
  [[ $a == out=* || $a == N=* ]] && copy=0
done
if [[ $copy == 1 ]]; then
  for f in fActiveMeanAndStd.pdf f2.png f2_high.png; do
    [[ -s "$here/$f" ]] || { echo "missing $here/$f, nothing copied" >&2; exit 1; }
  done
  for f in fActiveMeanAndStd.pdf f2.png f2_high.png; do
    cp "$here/$f" "$here/../../manuscript/figures/$f"
  done
  echo "copied fActiveMeanAndStd.pdf, f2.png, f2_high.png to manuscript/figures/"
fi
