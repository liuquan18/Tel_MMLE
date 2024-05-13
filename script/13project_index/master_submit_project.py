#!/usr/bin/env python
#%%
import os

#%%
models = ['MPI_GE','MPI_GE_onepct', 'CanESM2','CESM1_CAM5', 'GFDL_CM3','MK36']

#%%
for n, model in enumerate(models): #
    os.system(f"sbatch ./submit_project.sh {n+1} {model}")

# %%
