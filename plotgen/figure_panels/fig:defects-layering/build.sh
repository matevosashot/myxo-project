#!/usr/bin/env bash
# Builds defects-layering.pdf (fig:defects-layering) from its panels and copies it to
# manuscript/figures. Panels: (a) plotgen/defect_illustration/defect_{plus,minus}.pdf,
# (b) the 2x2 video-frame SVG below, made in myxo-video/ by
#   python figure_svg.py crop300-434 --anchor 08m08s --cols 2 --time-size 92 \
#     --label-size 76 --stroke 6
# from the frames of "Layers with defects 2.mp4" (the defect markers are in the movie),
# (c) plotgen/layering_profile/layering_profile.pdf. Only assembles; the panels are not
# regenerated. The SVG is converted with rsvg-convert if installed, otherwise with
# headless Chrome (CHROME overrides its path). Intermediates go to .tmp/.
# Runs from any directory.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"
name=defects-layering
svg="$(cd ../../../myxo-video/selected/annotated/Layers_with_defects_2_crop300-434-184x202_0.5x_up4-3 && pwd)/crop300-434_4frames_2x2_06m35s-08m08s-09m41s-11m13s.svg"
chrome="${CHROME:-/Applications/Google Chrome.app/Contents/MacOS/Google Chrome}"

mkdir -p .tmp
if command -v rsvg-convert > /dev/null; then
  rsvg-convert -f pdf -o .tmp/frames.pdf "$svg"
else
  read -r w h < <(sed -n '/<svg /{s/.*width="\([0-9.]*\)" height="\([0-9.]*\)".*/\1 \2/p;q;}' "$svg")
  printf '<!doctype html><style>@page{size:%spx %spx;margin:0}html,body{margin:0}img{display:block;width:%spx;height:%spx}</style><img src="file://%s">\n' \
    "$w" "$h" "$w" "$h" "$svg" > .tmp/frames.html
  "$chrome" --headless --disable-gpu --no-pdf-header-footer --allow-file-access-from-files \
    --user-data-dir="$here/.tmp/chrome" --print-to-pdf="$here/.tmp/frames.pdf" \
    "file://$here/.tmp/frames.html" 2> /dev/null
  rm -rf .tmp/chrome
fi

pdflatex -interaction=nonstopmode -halt-on-error -output-directory=.tmp "$name.tex" \
  > /dev/null || { tail -n 30 ".tmp/$name.log"; exit 1; }
grep '^panelrow:' ".tmp/$name.log"
mv ".tmp/$name.pdf" .
cp "$name.pdf" ../../../manuscript/figures/
echo "built $name.pdf, copied to manuscript/figures/$name.pdf"
