# %%
import xarray as xr
import numpy as np
import pandas as pd
import mpi4py.MPI as MPI
import sys
import os
import glob


# %%
import src.blocking.block_index as block_index
#%%
import importlib
importlib.reload(block_index)

# %%
# === mpi4py ===
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()  #0,1,2,3,4,5,6,7,8,9
    npro = comm.Get_size()  #10
except:
    print('::: Warning: Proceeding without mpi4py! :::')
    rank = 0
    npro = 1

#%%
num = int(sys.argv[1])
f1 = int(sys.argv[2])
f2 = int(sys.argv[3])

#%%
months = ["Jun", "Jul", "Aug","ano_Jun", "ano_Jul", "ano_Aug"]
month = months[num-1]
odir = f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/zg_{month}/"
files = glob.glob(odir + "*.nc")
files.sort()

#%%
list_all_pros = [0]*npro
for nn in range(npro):
  list_all_pros[nn] = files[nn::npro]
steps = list_all_pros[rank]



def detect_index(file):
        # open the file
        ds = xr.open_dataset(file)
        Z = ds.geopoth.squeeze()
        blocks = block_index.IB_index(Z,LLB_filter=True)
        return blocks

#%%
for kk, step in enumerate(steps):
    fname = os.path.basename(step)
    print(f'{month} on node: {num}: kk = {kk+1}/{len(steps)}, step = {fname}')
    to_path = f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_{month}/{fname}"
    blocks = detect_index(step)
    blocks.to_netcdf(to_path)

# %%
