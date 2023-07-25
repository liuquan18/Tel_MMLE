#!/usr/bin/env python
#%%
import os
import numpy as np

#%%
nsteps = 120
npar   = 12
njobs  = int(nsteps/npar) #10

#%%
for kk in range(njobs): #0,1,2
  k1 = kk*npar #0,4,8
  k2 = (kk+1)*npar #4,8,12
  os.system(f"sbatch ./submit_python_mpi4py.sh {kk+1} {k1} {k2}") # 1,0,4, 2,4,8, 3,8,12
# %%
