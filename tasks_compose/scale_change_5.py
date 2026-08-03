#!/usr/bin/env python3


import sys
import time
from pathlib import Path
import shutil
import os
import numpy as np

from myxo import FourierSolver
from myxo.correlation_def_v3 import GcalcForQ
import toolbox




SCALES = np.logspace(np.log10(100.0), np.log10(0.01), 200+1)
np.random.default_rng(801).shuffle(SCALES)
print(SCALES)


DATA_PATH = Path("/data/biophys/ashmat/data/myxo-project/scale_change_5.2/")
SCRATCH_PATH = os.environ.get("SCRATCH_PATH", "./")
SCRATCH_PATH = "/scratch/ashot/"

JOB_ID = os.environ.get("JOB_ID", None)
TASK_ID = os.environ.get("SLURM_ARRAY_TASK_ID", None)
if TASK_ID is None:
    print("No task ID provided. setting task ID to 0")
    TASK_ID = 0
else:
    TASK_ID = int(TASK_ID)

os.makedirs(f"/home/ashmat/Projects/myxo-project/tasks/logs/", exist_ok=True)
toolbox.setup_loggers(
        base_path=f"/home/ashmat/Projects/myxo-project/tasks/logs/{JOB_ID}.log",
        debug=True,
        stdout=True,
        train_logger=False,
    )

scale = SCALES[TASK_ID]

os.makedirs(DATA_PATH, exist_ok=True)

print(f"Running for scale {scale}")

model_params = {
    "S0": 1.0,
    "l": 1.0,
    "lM": 7.0 * scale,
    "lm": 0.5 * scale,
    "R": 16.0,
    "sigma": 6.0,
}

solver = FourierSolver(n=601, L=60, model_params=model_params,
                        dtype="float32", verbose=True, drop_phi=True, gcalc_Q_cls=GcalcForQ)
solver.precompute()

out = solver.save_results_h5(SCRATCH_PATH, save_P_diag=True, save_Q_diag=True)


print(f"Saved {out}")
_move_start = time.perf_counter()
shutil.move(out, DATA_PATH)
_move_elapsed = time.perf_counter() - _move_start
print(f"Moved {out} to {DATA_PATH} in {_move_elapsed:.2f}s")

