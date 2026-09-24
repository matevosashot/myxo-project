import numpy as np
from pack import *
rng=np.random.default_rng(0)
q=np.column_stack([rng.uniform(-.5,.5,(4,2)),rng.uniform(-3,3,4),rng.uniform(-3,3,4),rng.uniform(-10,10,4)])
q[:,:2]*=0.4
ell=np.full(4,ELL); fixed=np.zeros(4,bool)
E,g=energy(q,ell,fixed)
num=np.zeros_like(q); eps=1e-6
for i in range(4):
  for j in range(5):
    qp=q.copy();qp[i,j]+=eps;qm=q.copy();qm[i,j]-=eps
    num[i,j]=(energy(qp,ell,fixed,want_grad=False)-energy(qm,ell,fixed,want_grad=False))/2/eps
print(E); print(np.abs(num-g).max(), np.abs(g).max())
