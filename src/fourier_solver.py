from .correlation_def import Gcalc
import numpy as np
from scipy.fft import rfftn, irfftn, set_workers
import os




def _build_wavevectors_1d(n: int, dx: float) -> np.ndarray:
    return np.fft.fftfreq(n, d=dx / (2.0 * np.pi))


def _build_rfft_wavevectors_1d(n: int, dx: float) -> np.ndarray:
    return np.fft.rfftfreq(n, d=dx / (2.0 * np.pi))



def apply_fft(x_rank4):
    fft_axes = (0, 1, 2, 3)
    with set_workers(os.cpu_count()):
        comp_hat = rfftn(x_rank4, axes=fft_axes)
    return comp_hat

def apply_ifft(q_rank4):
    fft_axes = (0, 1, 2, 3)
    with set_workers(os.cpu_count()):
        x_rank4 = irfftn(q_rank4, axes=fft_axes)
    return x_rank4
    


class FourierSolver:
    def __init__(self, n, L,  model_params, dtype="float32", eps=None):
        self.n = n
        self.L = L
        self.model_params = model_params
        self.dtype = np.dtype(dtype)

        self.dx = 2.0 * L / n
        if eps is None:
            self.eps = (np.pi / (self.n * self.L)) ** 4
        else:
            self.eps = eps

        self.grid_1d = np.linspace(-L, L - self.dx, n)
        
        self.gcalc = Gcalc(self.grid_1d, **self.model_params, dtype=self.dtype)

    def precompute(self):
        self.gcalc.precompute()

        self.qfull = _build_wavevectors_1d(self.n, self.dx).astype(self.dtype)
        self.qhalf = _build_rfft_wavevectors_1d(self.n, self.dx).astype(self.dtype)

        self.qx  = self.qfull[:, None, None, None]
        self.qy  = self.qfull[None, :, None, None]
        self.qxp = self.qfull[None, None, :, None]
        self.qyp = self.qhalf[None, None, None, :]

        self.q1_inv_sq = 1.0 / (self.qx ** 2 + self.qy ** 2 + self.eps)
        self.q2_inv_sq = 1.0 / (self.qxp ** 2 + self.qyp ** 2 + self.eps)

        self.component = {0: self.qx, 1: self.qy}
        self.component_prime = {0: self.qxp, 1: self.qyp}

    def get_C_P_fft(self, a, b):
        c = self.gcalc.calc_C_P_index(a, b)

    def get_C_Q_fft(self, a, b, m, n):
        c = self.gcalc.calc_C_Q_index(a, b, m, n)
        

    def _k2_kernel(self, a, b):
        k1, k2 = self.component[a], self.component_prime[b]
        return k1 * k2 * self.q1_inv_sq * self.q2_inv_sq

    def _k4_kernel(self, a, b, m, n):
        k1, k2 = self.component[a], self.component_prime[b]
        k3, k4 = self.component[m], self.component_prime[n]
        return k1 * k2 * k3 * k4 * self.q1_inv_sq * self.q2_inv_sq

    def solve_P_contribution(self):
        s_pressure = None
        for a in range(2):
            for b in range(2):
                c = self.gcalc.calc_C_P_index(a, b)
                s = apply_fft(c) * self._k2_kernel(a, b)
                if s_pressure is None:
                    s_pressure = s
                else:
                    s_pressure += s
        s_pressure[0, 0, 0, 0] = 0.0
        
        c_pressure = apply_ifft(s_pressure)

        return c_pressure

    def solve_Q_contribution(self):
        s_pressure = None
        for a in range(2):
            for b in range(2):
                for m in range(2):
                    for n in range(2):
                        c = self.gcalc.calc_C_Q_index(a, b, m, n)
                        s = apply_fft(c) * self._k4_kernel(a, b, m, n)
                        if s_pressure is None:
                            s_pressure = s
                        else:
                            s_pressure += s
        c_pressure = apply_ifft(s_pressure)

        return c_pressure

        