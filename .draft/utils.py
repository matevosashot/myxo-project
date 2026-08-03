
import numpy as np
from concurrent.futures import ThreadPoolExecutor


def threaded_prod(A, B, out=None, workers=0):
    """
    Returns A * B using multiple threads. 
    """

    if out is None:
        out = np.empty((A.shape[0], B.shape[1]), dtype=A.dtype)

    def worker(idx):
        out[idx] = A[idx] @ B

    with ThreadPoolExecutor(max_workers=workers) as executor:
        executor.map(worker, range(A.shape[0]))

    return out