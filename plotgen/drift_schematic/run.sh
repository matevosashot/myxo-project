#!/bin/bash
# builds drift_schematic.pdf (illustration from slides4, frames m-drift and m-Y) and copies it to the manuscript
cd "$(dirname "$0")"
mkdir -p .tmp
pdflatex -interaction=nonstopmode -output-directory=.tmp drift_schematic.tex > /dev/null
cp .tmp/drift_schematic.pdf drift_schematic.pdf
cp drift_schematic.pdf ../../manuscript/figures/drift_schematic.pdf
