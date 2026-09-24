#!/bin/bash
cd "$(dirname "$0")"
PY=../.venv/bin/python
nohup $PY run_fill2.py brick3.pkl d1 10 0.3 300 > log_d1.txt 2>&1 &
nohup $PY run_fill2.py brick1.pkl d2 10 0.3 300 > log_d2.txt 2>&1 &
nohup $PY run_fill2.py brick3.pkl d3 3 0.3 300 > log_d3.txt 2>&1 &
nohup $PY run_fill2.py brick0.pkl d4 10 1.0 300 > log_d4.txt 2>&1 &
