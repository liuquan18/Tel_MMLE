import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import mpi4py.MPI as MPI
import glob
import os
import sys

# %%
def events_count(arr, threshold=10):
    """
    count the number of events that lasted more than 10 days for one single pixel
    """
    df = arr.to_dataframe()
    df = df['IB index'].reset_index()
    G=df[df['IB index']> 0].groupby((df['IB index']<=0).cumsum())
    durations = G.size()
    count = durations[durations> threshold].count()
    count_x = xr.DataArray(count, dims = ['z'],coords = {'z':arr.z})
    return count_x

def yearly_events_count(arr, threshold=10):
    """
    statistcis (average, count) of the duration of events 
    """
    arr = arr.stack(z=('lat', 'lon'))
    count = arr.groupby('z').apply(events_count, threshold=threshold)
    return count.unstack()

def decade_events_count(arr, threshold=10):
    decade_count = arr.resample(time = '10AS').apply(yearly_events_count, threshold=threshold)
    return decade_count


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

anomaly = False

#%%
if anomaly:
    odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_event_ano/"
    todir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/block_event_count_ano/"
else:
    odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_event/"
    todir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/block_event_count/"

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
