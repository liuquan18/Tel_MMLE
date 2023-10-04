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
import src.blocking.wave_breaking as wave_breaking

# %%
import importlib

importlib.reload(block_index)

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

anomaly = True

#%%
var_name = "wb"


def detect_index(file, var_name="block"):
    # open the file
    ds = xr.open_dataset(file)
    Z = ds.geopoth.squeeze()
    if var_name == "block":
        blocks = block_index.IB_index(Z, LLB_filter=False)
    elif var_name == "wb":
        blocks = wave_breaking.wave_breaking_index(Z)
    return blocks


# %%
if anomaly:
    month = ["ano_Jun", "ano_Jul", "ano_Aug"]
else:
    month = ["Jun", "Jul", "Aug"]
files_globes = []
for mm in month:
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/zg_{mm}/"
    files_globe = glob.glob(odir + "*.nc")
    files_globe.sort()
    files_globes.append(files_globe)
# flat the files_globes into a list
files_globes = [item for sublist in files_globes for item in sublist]
files = files_globes[f1:f2]

# %%
list_all_pros = [0] * npro
for nn in range(npro):
    list_all_pros[nn] = files[nn::npro]
steps = list_all_pros[rank]


# %%
for kk, step in enumerate(steps):
    fname = os.path.basename(step)
    print(f"{month} on node: {num}: kk = {kk+1}/{len(steps)}, step = {fname}")
    # replace the 'zg' in step with 'block'
    to_path = step.replace("zg_", f"{var_name}_")
    blocks = detect_index(step, var_name=var_name)
    blocks.to_netcdf(to_path)

# %%
