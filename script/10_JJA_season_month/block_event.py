
#%%
import xarray as xr
import numpy as np
import src.blocking.block_event as block_event
from src.blocking.utils import reg_lens
import mpi4py as MPI
import glob
import os
import sys

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
#%%
odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_"

# Create an array of integers from 69 to 100 (inclusive)
all_ens = np.arange(69, 101, 1)

# Convert the integers to strings with leading zeros
all_ens_str = [str(i).zfill(4) for i in all_ens]

global_files = np.empty((32,3), dtype=object)
for i, ens in enumerate(all_ens_str):
    fname_pre = f"onepct_1850-1999_ens_{ens}.gph500"
    if anomaly:
        fnames = [odir + 'ano_Jun/' + fname_pre + '_06.nc', 
                    odir + 'ano_Jul/' + fname_pre + '_07.nc', 
                    odir + 'ano_Aug/' + fname_pre + '_08.nc']
    else:
        fnames = [odir + 'Jun/' + fname_pre + '_06.nc', 
            odir + 'Jul/' + fname_pre + '_07.nc', 
            odir + 'Aug/' + fname_pre + '_08.nc']
    global_files[i,:] = fnames

# %%
files = global_files[f1:f2]
list_all_pros = [0] * npro
for nn in range(npro):
    list_all_pros[nn] = files[nn::npro]
steps = list_all_pros[rank]

#%%
if anomaly:
    todir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_event_ano/"
else:
    todir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_event/"

for kk, step in enumerate(steps):
    print(f"node {num} Process {rank} is working on {kk+1}/{len(steps)}")
    inds = []
    for i, file in enumerate(step):
        ds = xr.open_dataset(file)
        inds.append(ds['IB index'])
    ix = xr.concat(inds, dim = 'time').sortby('time')

    ix_box = block_event._box_ix(ix)
    BE=block_event._yearly_persistence(ix_box,5)
    BE.to_netcdf(todir + step[0][-38:-6] + '.nc')
    print(f"Process {rank} is done with {kk+1}/{len(steps)}")
# %%
