#!/bin/bash
cd "$(dirname "$0")"
PY=../.venv/bin/python
nohup $PY run_flex.py state_d1.pkl f1 0.5 > log_f1.txt 2>&1 &
nohup $PY run_flex.py state_d4.pkl f4 0.5 > log_f4.txt 2>&1 &
nohup $PY run_flex.py state_d1.pkl f1b 0.05 > log_f1b.txt 2>&1 &
nohup $PY run_flex.py state_d4.pkl f4b 0.05 > log_f4b.txt 2>&1 &
