#!/bin/bash
cd "$(dirname "$0")"
PY=../.venv/bin/python
nohup $PY run_fill.py b10_mid.pkl c1 10 0.5 0.3 > log_c1.txt 2>&1 &
nohup $PY run_fill.py b10_mid.pkl c2 10 0.5 1.0 > log_c2.txt 2>&1 &
nohup $PY run_fill.py seed0.pkl c3 10 0.5 0.3 > log_c3.txt 2>&1 &
