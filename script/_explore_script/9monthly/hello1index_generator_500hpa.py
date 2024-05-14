# %%
# %%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as np
import sys
import time as pytime
import pandas as pd

# %%
import importlib

importlib.reload(index_generate)

#%%
num = int(sys.argv[1])
t1 = int(sys.argv[2])
t2 = int(sys.argv[3])
#%%
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
models = ["GFDL"]#["MPI_GE","CanESM2","CESM1_CAM5","GFDL_CM3","MK36"] #,"GFDL_CM3"
# the steps that need to be run on this node
steps_global = np.arange(len(months))
steps = steps_global[t1:t2]
model = models[num-1]
# model = "CanESM2"


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


# %%
# function for generate the index
def index_gen(season,model=model):
    generator = index_generate.decompose_plev(
        model =model, plev=50000, fixedPattern='decade', standard='first', season=season
    )
    generator.save_result()

#%%
# use mpi4py
for kk, step in enumerate(steps):
    print("**========**")
    print(f'node: {num}: kk = {kk+1}/{steps.size}, model = {model}, month = {months[step]}')
    index_gen(months[step])

# %%
