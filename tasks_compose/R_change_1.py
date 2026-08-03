#!/usr/bin/env python3


import sys
import toolbox
from pathlib import Path
import shutil
import os
import numpy as np

from myxo import FourierSolver
from myxo.correlation_def_v2_2 import GcalcForQ



Rs = np.linspace(1,80, 41)
np.random.default_rng(801).shuffle(Rs)
print(Rs)


DATA_PATH = Path("/data/biophys/ashmat/data/myxo-project/R_change_1/")
SCRATCH_PATH = os.environ.get("SCRATCH_PATH", "./")
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

R = Rs[TASK_ID]

os.makedirs(DATA_PATH, exist_ok=True)

print(f"Running for R {R}")

model_params = {
    "S0": 1.0,
    "l": 1.0,
    "lM": 7.0,
    "lm": 0.5,
    "R": R,
    "sigma": 6.0,
}

solver = FourierSolver(n=601, L=60, model_params=model_params,
                        dtype="float32", verbose=True, drop_phi=True, gcalc_Q_cls=GcalcForQ)
solver.precompute()

out = solver.save_results_h5(SCRATCH_PATH, save_P_diag=True, save_Q_diag=True)
print(f"Saved {out}")
shutil.move(out, DATA_PATH)
print(f"Moved {out} to {DATA_PATH}")

