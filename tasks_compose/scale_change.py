#!/usr/bin/env python3


import sys
from pathlib import Path
import shutil
import os
import numpy as np

from myxo import FourierSolver
import toolbox




SCALES = np.linspace(0.2, 2.0, 41)[::-1]


DATA_PATH = Path("/data/biophys/ashmat/data/myxo-project/scale_change_2/")
SCRATCH_PATH = os.environ["SCRATCH_PATH"]
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

solver = FourierSolver(n=500, L=50, model_params=model_params,
                        dtype="float32", verbose=True)
solver.precompute()

out = solver.save_results_h5(SCRATCH_PATH, save_P_diag=True, save_Q_diag=True)
print(f"Saved {out}")
shutil.move(out, DATA_PATH)
print(f"Moved {out} to {DATA_PATH}")

