import os

_NWORKERS = len(os.sched_getaffinity(0))
os.environ["NUMEXPR_MAX_THREADS"] = str(_NWORKERS)

import numexpr as ne
ne.set_num_threads(_NWORKERS)

import toolbox
toolbox.setup_loggers(
    base_path="./myxo.log",
    debug=False,
    stdout=True,
    train_logger=False,
)



from .correlation_def import Gcalc
from .sparse_tensor import SparseTensor
from .fourier_solver import FourierSolver