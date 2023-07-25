# %%
import src.MMLE_TEL.story_line as story_line
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf
import src.MMLE_TEL.index_stats as index_stats
import xarray as xr
import os
import sys
#%%
import importlib

importlib.reload(story_line)
importlib.reload(extrc_tsurf)


#%%
num = int(sys.argv[1])
t1 = int(sys.argv[2])
t2 = int(sys.argv[3])
#%%
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
models = ["MPI_GE"]#["MPI_GE","CanESM2","CESM1_CAM5","MK36","GFDL_CM3"]
# the steps that need to be run on this node
steps_global = np.arange(len(months))
steps = steps_global[t1:t2]
model = models[num-1]

#%%
# === mpi4py ===
try:
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  npro = comm.Get_size()
except:
  print('::: Warning: Proceeding without mpi4py! :::')
  rank = 0
  npro = 1


list_all_pros = [0]*npro
for nn in range(npro):
  list_all_pros[nn] = steps[nn::npro]
steps = list_all_pros[rank]


#%%
# MPI_GE_onepct
# MPI_GE_onepct summer
def extreme_count(season,model = model):
    index_stats.extreme_counts_tsurf(
        model        = model,
        tsurf        =  "ens_fld_year_mean", 
        standard     =  "first",
        season       =  season,)
    

#%%
# use mpi4py
for kk, step in enumerate(steps):
    print("**========**")
    print(f'node: {num}: model = {model}, kk = {kk+1}/{steps.size}, month = {months[step]}')
    extreme_count(months[step])


