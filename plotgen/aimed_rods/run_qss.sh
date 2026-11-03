#!/usr/bin/env bash
# Regenerates aimed_rods_qss.pdf and aimed_rods_qss.png (fig:aimed-rods of the manuscript,
# Sec. IV since 2026-10-23: rods aimed at, or side-on to, a point in a uniform nematic, at
# the extremes of the variance on the head and tail sides and at the core of the moving +1/2 defect, on its
# steady state Q_ss as the dotted curves of fig:variance) next to this script, and copies the
# pdf to manuscript/figures/. Four panels (uniform, head side, core, tail side) for \textwidth.
# Arguments are passed on (nocore, ssDir, Rc=<um>; out=<dir> writes elsewhere and skips the
# copy; nocore writes the earlier three panels, aimed_rods_qss_3panel.pdf, and copies nothing). Runs from any directory.
set -eo pipefail                   # no -u: lmod.sh reads unset variables

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"

type module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh   # non-login shells
module load mathematica/14.10
unset DISPLAY                       # headless export
export QT_QPA_PLATFORM=offscreen

wolframscript -file "$here/aimed_rods_qss.wls" "$@"

if [[ "$*" != *out=* && " $* " != *" nocore "* ]]; then   # nocore: the old 3 panels, not copied
  cp "$here/aimed_rods_qss.pdf" "$repo/manuscript/figures/aimed_rods_qss.pdf"
  echo "copied aimed_rods_qss.pdf to manuscript/figures/"
fi
