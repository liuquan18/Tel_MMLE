#!/usr/bin/env python
#%%
import os
import numpy as np

#%%
nsteps = 30
npar   = 2
njobs  = int(nsteps/npar) # 15

#%%
for kk in range(njobs): #
    k1 = kk*npar #
    k2 = (kk+1)*npar #
    os.system(f"sbatch ./block_generator_submitter.sh {kk+1} {k1} {k2}") # 1,0,4, 2,4,8, 3,8,12

# %%
