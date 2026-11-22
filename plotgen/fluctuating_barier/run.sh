#!/usr/bin/env bash
# Builds fluctuating_barrier.pdf (fig:nucleation(b)) next to this script. Runs from any
# directory; extra arguments go to the Python script (e.g. --n 200 --seed 3).
# plotgen/figure_panels/fig:nucleation/build.sh places it in the row; run that next.
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py="$here/../../.venv/bin/python"   # system python3 has no matplotlib

"$py" "$here/plot_fluctuating_barier.py" --n 500 "$@"
