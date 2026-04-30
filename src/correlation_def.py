import argparse
from termios import PARODD
import numpy as np
import numexpr as ne

def f(x):
    """f[x_] := x Sqrt[(0.34 + 0.07 x^2)/(1 + 0.41 x^2 + 0.07 x^4)]"""
    x2 = x**2
    return x * np.sqrt((0.34 + 0.07 * x2) / (1.0 + 0.41 * x2 + 0.07 * x2**2))

def S(r, S0, l, R):
    """
    Strength of the nematic order parameter.
    S[r_] := S0 f[r/l]
    """
    return S0 * f(r / l) * np.exp(-(r / R))  # add Gaussian cutoff to ensure finite integrals

def _mvn_pdf_2d(dx, dy, cov):
    """
    Vectorised 2-D multivariate normal PDF from displacement components.
    dx, dy : broadcastable scalar arrays  (displacement r2 - r1)
    cov    : (..., 2, 2)  — may have fewer broadcast dims than dx/dy
    returns: broadcast shape of (dx, dy, cov[..., 0, 0])
    """
    dtype = cov.dtype
    det  = cov[..., 0, 0] * cov[..., 1, 1] - cov[..., 0, 1] * cov[..., 1, 0]
    # Mahalanobis: uses only 3 unique inverse entries (cov is symmetric)
    maha2 = (  cov[..., 1, 1] / det * dx**2
             - dtype.type(2) * cov[..., 0, 1] / det * dx * dy
             +  cov[..., 0, 0] / det * dy**2 )
    return np.exp(dtype.type(-0.5) * maha2) / (dtype.type(2 * np.pi) * np.sqrt(det))

def _mvn_pdf_2d_fast(dx, dy, cov):
    c00 = cov[..., 0, 0]
    c01 = cov[..., 0, 1]
    c11 = cov[..., 1, 1]
    t = cov.dtype.type
    half = t(-0.5)
    two = t(2)
    inv_two_pi = t(1.0 / (2 * np.pi))
    return ne.evaluate(
        "exp(half * (c11 * dx**2 - two * c01 * dx * dy + c00 * dy**2) / (c00 * c11 - c01 * c01))"
        " * inv_two_pi / sqrt(c00 * c11 - c01 * c01)"
    )




class Gcalc:
    def __init__(self, grid_1d, S0=1.0, l=1.0, lM=7.0/3, lm=0.7/3, R=10, dtype=np.float64):
        self.S0 = S0
        self.l = l
        self.lM = lM
        self.lm = lm
        self.R = R

        self.grid = grid_1d.astype(dtype)
        self.x1 = self.grid[:, None, None, None]
        self.y1 = self.grid[None, :, None, None]
        self.x2 = self.grid[None, None, :, None]
        self.y2 = self.grid[None, None, None, :]

        dx = grid_1d[1] - grid_1d[0]
        self.L_grid = grid_1d[-1] - grid_1d[0] + dx  # total length of the periodic box
    
        self.dtype = dtype


    def calcQ(self):
        """
        Q tensor (2x2) at positions (x, y).
        x, y can be arrays of shape (...).
        Returns shape (..., 2, 2).
        """
        x = self.grid[:, None]
        y = self.grid[None, :]

        dtype = np.result_type(x, y)
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)          # Mathematica ArcTan[x, y] = atan2(y, x)
        theta = phi / dtype.type(2)     # Theta[phi] := phi/2
        s = S(r, self.S0, self.l, self.R)                 # scalar field, shape (...)
        c2 = np.cos(dtype.type(2) * theta)
        s2 = np.sin(dtype.type(2) * theta)

        # Build (..., 2, 2) tensor:  [[S cos2θ, S sin2θ], [S sin2θ, -S cos2θ]]
        out = np.empty(np.broadcast_shapes(np.shape(x), np.shape(y)) + (2, 2), dtype=dtype)
        out[..., 0, 0] = s * c2
        out[..., 0, 1] = s * s2
        out[..., 1, 0] = s * s2
        out[..., 1, 1] = -s * c2

        self.Q = out

    def calcP(self):
        """P = 1/2 (I + Q)  — shape (..., 2, 2)"""
        Q = self.Q
        eye = np.eye(2, dtype=self.dtype)
        self.P = self.dtype.type(0.5) * (eye + Q)

    def calcSigma(self):
        """
        Σ = 1/2 I (lM² + lm²) + 1/2 Q (lM² - lm²)
        Returns shape (..., 2, 2) covariance matrix.
        """
        dtype = self.dtype
        Q = self.Q
        eye = np.eye(2, dtype=dtype)
        self.Sigma = 0.5 * eye * (self.lM**2 + self.lm**2) + 0.5 * Q * (self.lM**2 - self.lm**2)
   
    def calc_phi(self):
        Q = self.Q 
        Sigma = self.Sigma

        dx = self.x2 - self.x1
        dx = dx - self.L_grid * np.round(dx / self.L_grid)  # periodic boundary conditions

        dy = self.y2 - self.y1
        dy = dy - self.L_grid * np.round(dy / self.L_grid)  # periodic boundary conditions

        self.phi = _mvn_pdf_2d_fast(dx, dy, Sigma[:, :, None, None])

        self.phi_transpose = self.phi.transpose(2, 3, 0, 1)

    def calc_C_P_index(self, a, b):
        """
        Calculate the C_Q tensor (2x2) at positions (x, y).
        x, y can be arrays of shape (...).
        Returns shape (..., 2, 2).
        """
        phi = self.phi # N x N x N x N
        phi_T = self.phi.transpose(2, 3, 0, 1)
        P = self.P[:, :, None, None] # N x N x 1 x 1 x 2 x 2
        PT = P.transpose(2, 3, 0, 1, 4, 5) # 1 x 1 x N x N x 2 x 2

        return ne.evaluate(
            "half * ( phi * P_ab + phi_T * PT_ab)",
            local_dict={
                'half': self.dtype.type(0.5),
                'phi': phi, 'phi_T': phi_T, 
                'P_ab': P[..., a, b], 'PT_ab': PT[..., a, b]
            })

    def calc_C_P(self):
        """
        Calculate the C_P tensor (2x2) at positions (x, y).
        x, y can be arrays of shape (...).
        Returns shape (..., 2, 2).
        """
        N = len(self.grid)

        values = ne.evaluate("half * (P_left * phi_left + P_right * phi_right)", 
            local_dict={
                'half': self.dtype.type(0.5), 
                'P_left': self.P[:,:,None,None], 'phi_left': self.phi[..., None, None], 
                'P_right': self.P[None, None, :, :], 'phi_right': self.phi_transpose[..., None, None]
            })

        # else:
        #     values = ne.evaluate("half * (P_left * phi_left + P_right * phi_right)" , 
        #         local_dict={
        #             'half': self.dtype.type(0.5), 
        #             'P_left': self.P[:,:,None,None, i, j], 'phi_left': self.phi,
        #             'P_right': self.P[None, None, :, :, i, j], 'phi_right': self.phi_transpose
        #         })

        return values


    def calc_C_Q_index(self, a, b, m, n):
        phi = self.phi                                    # N x N x N x N
        phi_T = self.phi_transpose                        # N x N x N x N
        P = self.P[:, :, None, None]                      # N x N x 1 x 1 x 2 x 2
        PT = P.transpose(2, 3, 0, 1, 4, 5)                # 1 x 1 x N x N x 2 x 2

        return ne.evaluate(
            "half * phi   * (P_ab  * P_mn  + P_am  * P_bn  + P_an  * P_bm  - P_ab * PT_mn)"
            " + "
            "half * phi_T * (PT_ab * PT_mn + PT_am * PT_bn + PT_an * PT_bm - P_ab * PT_mn)",
            local_dict={
                'half': self.dtype.type(0.5),
                'phi': phi, 'phi_T': phi_T,
                'P_ab': P[..., a, b], 'P_mn': P[..., m, n],
                'P_am': P[..., a, m], 'P_bn': P[..., b, n],
                'P_an': P[..., a, n], 'P_bm': P[..., b, m],
                'PT_ab': PT[..., a, b], 'PT_mn': PT[..., m, n],
                'PT_am': PT[..., a, m], 'PT_bn': PT[..., b, n],
                'PT_an': PT[..., a, n], 'PT_bm': PT[..., b, m],
            })

    def calc_C_Q(self, a=None, b=None, m=None, n=None):
        """
        Calculate the C_Q 2 dimensional rank 4 tensor (2x2x2x2) at positions (x, y) on the grid.
        Returns shape (N, N, N, N, 2, 2, 2, 2).
        """
        N = len(self.grid)
        
        phi = self.phi # N x N x N x N
        P = self.P[:, :, None, None] # N x N x 1 x 1 x 2 x 2
        P_transpose = P.transpose(2,3,0,1,4,5) # 1 x 1 x N x N x 2 x 2

        # i in the j-th place in the tuple means that the array's i-th axis becomes the transposed array's j-th axis.
        #            4567 
        # P_outer_rr_abmn = P[..., :, :, None, None] * P[..., None, None, :, :]
        P_outer_rr_abmn = ne.evaluate("P_left * P_right", 
            local_dict={
                'P_left': P[..., :, :, None, None], 
                'P_right': P[..., None, None, :, :]
            }) # N x N x 1 x 1 x 2 x 2 x 2 x 2
        P_outer_rr_ambn = P_outer_rr_abmn.transpose(0,1,2,3,   4, 6, 5, 7)
        P_outer_rr_anbm = P_outer_rr_abmn.transpose(0,1,2,3,   4, 6, 7, 5)

        # P_outer_rr1_abmn = P[...,:,:,None,None] * P_transpose[...,None,None,:,:] # N x N x N x N x 2 x 2 x 2 x 2
        P_outer_rr1_abmn = ne.evaluate("P_left * P_right", 
            local_dict={
                'P_left': P[...,:,:,None,None], 
                'P_right': P_transpose[...,None,None,:,:]
            }) # N x N x N x N x 2 x 2 x 2 x 2
        
        # part_left = phi[...,None,None,None,None] * (P_outer_rr_abmn + P_outer_rr_ambn + P_outer_rr_anbm - P_outer_rr1_abmn)
        part_left = ne.evaluate("half * phi * (P_outer_rr_abmn + P_outer_rr_ambn + P_outer_rr_anbm - P_outer_rr1_abmn)", 
            local_dict={
                'half': self.dtype.type(0.5), 
                'phi': phi[...,None,None,None,None], 
                'P_outer_rr_abmn': P_outer_rr_abmn, 
                'P_outer_rr_ambn': P_outer_rr_ambn, 
                'P_outer_rr_anbm': P_outer_rr_anbm, 
                'P_outer_rr1_abmn': P_outer_rr1_abmn
            }) # N x N x N x N x 2 x 2 x 2 x 2
        #  0123 4567 
        #  2301 6745
        part_right = part_left.transpose(2,3,0,1,6,7,4,5)

        value = ne.evaluate("part_left + part_right", 
            local_dict={
                'part_left': part_left, 
                'part_right': part_right
            }) # N x N x N x N x 2 x 2 x 2 x 2



        return value

        
    def precompute(self):
        self.calcQ()
        self.calcP()
        self.calcSigma()
        self.calc_phi()

    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute G_{ij} on a 4D grid")
    parser.add_argument("--S0",       type=float, default=1.0)
    parser.add_argument("--l",        type=float, default=1.0)
    parser.add_argument("--lM",       type=float, default=7.0/3)
    parser.add_argument("--lm",       type=float, default=0.7/3)
    parser.add_argument("--R",        type=float, default=10, help="Cutoff radius")
    parser.add_argument("--N",        type=int,   default=50, help="Grid size per dimension")
    parser.add_argument("--parallel", action="store_true",     help="Compute G in parallel using threads")
    parser.add_argument("--workers",  type=int,   default=None, help="Number of worker threads (default: CPU count)")
    parser.add_argument("--dtype",    type=str,   default="float64", help="NumPy dtype for computation (e.g. float16, float32, float64)")
    args = parser.parse_args()

    dtype = np.dtype(args.dtype)

    print("O(N^4) memory is {:0.1f}GB".format(args.N**4 / 1024**3 * dtype.itemsize))



    grid_1d = np.linspace(-5, 5, args.N, dtype=dtype)
    gcalc = Gcalc(grid_1d, S0=args.S0, l=args.l, lM=args.lM, lm=args.lm, R=args.R, dtype=dtype)
    gcalc.precompute()
    print(gcalc.Sigma.mean())
    # print(gcalc.phi[4,10:12,3,5:7])
    # print(gcalc.calc_C_P(i=0, j=0).mean())
    # print(gcalc.calc_C_P().mean())
    value = gcalc.calc_C_P()

    value_2 = np.empty_like(value)
    for a in range(2):
        for b in range(2):
            value_2[..., a, b] = gcalc.calc_C_P_index(a, b)
    print("absdiff", np.abs(value - value_2).mean())
    

    value = gcalc.calc_C_Q()

    value_2 = np.empty_like(value)
    for a in range(2):
        for b in range(2):
            for m in range(2):
                for n in range(2):
                    value_2[..., a, b, m, n] = gcalc.calc_C_Q_index(a, b, m, n)

    print("absdiff", np.abs(value - value_2).mean())
    exit(0)