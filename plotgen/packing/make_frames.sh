#!/bin/bash


python render_packing.py --out frames/frame_1.pdf

python render_packing.py --out frames/frame_2.pdf \
    --show-streams

python render_packing.py --out frames/frame_3.pdf \
    --show-streams \
    --pop-cell
