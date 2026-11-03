#!/usr/bin/env bash
# Builds fluctuating_barrier.pdf (fig:nucleation(b)) next to this script. Runs from any
# directory; extra arguments go to the Python script (e.g. --n 200 --seed 3).
# The panel lives inside figures/handmade/main_v2.pdf, so it is not copied anywhere.
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py="$here/../../.venv/bin/python"   # system python3 has no matplotlib

"$py" "$here/plot_fluctuating_barier.py" --n 500 "$@"
