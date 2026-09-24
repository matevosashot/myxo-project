#!/bin/bash
cd "$(dirname "$0")"
PY=../.venv/bin/python
for sd in 0 1 2 3; do nohup $PY run_flex.py fseed$sd.pkl g$sd 0.5 > log_g$sd.txt 2>&1 & done
