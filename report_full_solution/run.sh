#!/bin/bash
# Rebuild the report from the solver output in ../comoving_defect_full/out/:
#   1. export_fields.wls   fields and diagnostics of every run -> data/*.json
#   2. dispersion.py       linear stability of the uniform state -> figures/dispersion.pdf
#   3. collect_numbers.py  numbers quoted in the text -> numbers.tex
#   4. make_figures.py     figures/*.pdf
#   5. fluctuations: export_fluct_ref.wls (the paper's ssDir runs -> data/fluct_ref_*.json),
#      fluct_analysis.py (../data/fluctuations-full/* -> data/fluct_summary.json),
#      fluct_numbers.py (numbers_fluct.tex), make_fluct_figures.py (figures/fluct_*.pdf)
#   6. latexmk             report.pdf
# conditioning.wls (data/conditioning.json) and check_jacobian.wls
# (data/check_jacobian.log) are run by hand; see report.tex Sec. 4-5.
set -eu
cd "$(dirname "$0")"
unset DISPLAY
WS=/usr/local/math/math1410/Executables/wolframscript
PY=../.venv/bin/python
TAGS="A1 B1 C1 Z1 Z15 H1 A15 F15 F2 G2 L10 L15 L30 L40"
avail=""
for t in $TAGS; do
  [ -s "../comoving_defect_full/out/$t.mx" ] && avail="$avail $t"
done
echo "runs found:$avail"
$WS -file export_fields.wls $avail
$PY -I dispersion.py
$PY -I collect_numbers.py > data/numbers.txt
$PY -I make_figures.py --A A1 --B B1
[ -s data/fluct_ref_N180.json ] || $WS -file export_fluct_ref.wls
$PY -I fluct_analysis.py > data/fluct_analysis.txt
$PY -I fluct_numbers.py > /dev/null
$PY -I make_fluct_figures.py
latexmk -pdf -interaction=nonstopmode -halt-on-error report.tex > data/latexmk.log 2>&1 \
  || { tail -n 40 data/latexmk.log; exit 1; }
echo "wrote report.pdf"
