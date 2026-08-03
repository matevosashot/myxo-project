import myxo
import numpy as np
from myxo import FourierSolver
from myxo import _NWORKERS

N = 400
L = 50
scale = 1.0
model_params = dict(
    S0=1.0,
    l=1.0,
    lM=7.0 * scale,
    lm=0.5 * scale,
    R=16,
    sigma=6.0,
)
solver = FourierSolver(n = N, L = L, model_params = model_params, verbose=True, mmap=False)

solver.precompute()

P_contrib = solver.solve_P_contribution__symmetry_optimized()


P_contrib_diag = P_contrib[np.arange(N)[:, None], np.arange(N)[None, :], 
               np.arange(N)[:, None], np.arange(N)[None, :]]




import plotly.graph_objects as go

fig = go.Figure(data=[go.Surface(z=P_contrib_diag, colorscale='viridis')])
fig.update_layout(scene=dict(xaxis_title='j', yaxis_title='i', zaxis_title='P_contrib_diag'))
fig.write_html("P_contrib_diag.html")
# fig.show()