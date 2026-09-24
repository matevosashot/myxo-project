import numpy as np, pickle
import flex
from flex import *
q,ell,f=pickle.load(open('state_d1.pkl','rb'))
V=from_q(q,ell)[:8]; ell=ell[:8]; fixed=np.zeros(8,bool)
rng=np.random.default_rng(1); V[:,2:]+=rng.normal(0,0.3,V[:,2:].shape); V[:,:2]+=0.0
pass
E,g=energy(V,ell,fixed); num=np.zeros_like(V); eps=1e-6
for i in range(8):
  for j in range(V.shape[1]):
    a=V.copy();a[i,j]+=eps;b=V.copy();b[i,j]-=eps
    num[i,j]=(energy(a,ell,fixed,want_grad=False)-energy(b,ell,fixed,want_grad=False))/2/eps
print(E, np.abs(num-g).max(), np.abs(g).max())
i,j=np.unravel_index(np.abs(num-g).argmax(),g.shape); print(i,j,num[i,j],g[i,j])
X,phi,h=geom(V,ell); print(np.linalg.norm(0.5*(X[i,:-1]+X[i,1:]),axis=1).min())
