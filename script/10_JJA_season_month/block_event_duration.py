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

# %%

def average_duration(event):
    """
    the average duraion for one single pixel
    """
    
    df = event.to_dataframe()
    df = df['IB index'].reset_index()
    # see https://mp.weixin.qq.com/s?__biz=MzkzNjI1ODQxMQ==&mid=2247484296&idx=1&sn=a9b63336bfbf13724d8b9d13321c75e4&chksm=c2a03d8cf5d7b49afef7ccbf4e0337d019ec1419a241a45ca8e6ca5cfcba6658dc39a9eebd04&token=1605307128&lang=zh_CN#rd
    Grouper = df.groupby(df.time.dt.year)['IB index'].transform(lambda x: (x<=0).cumsum())
    G=df[df['IB index']> 0].groupby([df.time.dt.year,Grouper])
    durations = G.size()
    dur_mean = durations.mean()
    dur_mean_x = xr.DataArray(dur_mean, dims = ['z'],coords = {'z':event.z})
    return dur_mean_x

def spatial_applyer(event):
    event = event.stack(z=('lat', 'lon'))
    mean_dur = event.groupby('z').apply(average_duration)


def decade_average_duration(arr):
    decade_dur = arr.resample(time = '10AS').apply(spatial_applyer)
    return decade_dur



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
nollb = False
#%%
o_pre = "block_nollb_event" if nollb else "block_event"
o_pre = o_pre + "_ano/" if anomaly else o_pre+"/"
#%%
to_pre = "block_nollb_duration_average" if nollb else "block_duration_average"
to_pre = to_pre + "_ano/" if anomaly else to_pre+"/"

#%%

odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/" + o_pre
todir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/" + to_pre

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
    average_dur = decade_average_duration(event)
    average_dur.name = 'average_duration'
    average_dur.to_netcdf(todir + 'dec_' + fbase_name)

# %%
