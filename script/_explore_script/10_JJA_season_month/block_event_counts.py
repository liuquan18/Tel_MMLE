#%%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import mpi4py.MPI as MPI
import glob
import os
import sys
from src.blocking.block_event_stat import decade_events_count

# %%
# === mpi4py ===
try:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()  # 0,1,2,3,4,5,6,7,8,9
    npro = comm.Get_size()  # 10
except:
    print("::: Warning: Proceeding without mpi4py! :::")
    rank = 0
    npro = 1

# %%
num = int(sys.argv[1])
f1 = int(sys.argv[2])
f2 = int(sys.argv[3])

#%%
anomaly = False
nollb = True
#%%
o_pre = "block_nollb_event" if nollb else "block_event"
o_pre = o_pre + "_ano/" if anomaly else o_pre+"/"
#%%
to_pre = "block_nollb_event_count" if nollb else "block_event_count"
to_pre = to_pre + "_ano/" if anomaly else to_pre+"/"

#%%

odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/" + o_pre
todir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/" + to_pre

#%%
# mkdir if todir does not exist
if not os.path.exists(todir):
    os.makedirs(todir)

global_files = glob.glob(odir + "*.nc")
files = global_files[f1:f2]
list_all_pros = [0] * npro
for nn in range(npro):
    list_all_pros[nn] = files[nn::npro]
steps = list_all_pros[rank]

#%%
for kk, step in enumerate(steps):
    fbase_name = os.path.basename(step)
    print(f"node {num} Process {rank} is working on {kk+1}/{len(steps)}")
    print(f"current file is {fbase_name}")
    event = xr.open_dataset(step)
    event_count = decade_events_count(event)
    event_count.name = 'block_event_count'
    event_count.to_netcdf(todir + 'dec_' + fbase_name)

# %%
