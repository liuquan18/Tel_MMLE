
#%%
import xarray as xr
import numpy as np
import src.blocking.block_event as block_event
from src.blocking.utils import reg_lens
import mpi4py as MPI
import glob
import os
import sys
import src.blocking.block_event_stat as block_event_stat

#%%
import importlib
importlib.reload(block_event)
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
prefix =  'block_nollb_event_' if nollb else 'block_event_'
folder = prefix + "ano_" if anomaly else prefix
#%%
odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/" + folder

months = ['Jun', 'Jul', 'Aug']
files_globes = []
for month in months:
    files_globe = glob.glob(odir + month + "/*.nc")
    files_globe.sort()
    files_globes.append(files_globe)

files_globes = [item for sublist in files_globes for item in sublist]
files = files_globes[f1:f2]

#%%
list_all_pros = [0] * npro
for nn in range(npro):
    list_all_pros[nn] = files[nn::npro]
steps = list_all_pros[rank]

# %%
for kk, step in enumerate(steps):
    fname = os.path.basename(step)
    print(f"{month} on node: {num}: kk = {kk+1}/{len(steps)}, step = {fname}")
    # replace the 'zg' in step with 'block'
    to_path = step.replace(prefix, prefix[:-6]+"average_duration_")
    to_path = to_path.replace("30_daily", "30")
    ix = xr.open_dataset(step)['IB index']
    BE = block_event_stat.annual_average_duration(ix) # 10 days (6h a day)
    BE.to_netcdf(to_path)

# %%
