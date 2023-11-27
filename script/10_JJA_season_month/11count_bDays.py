#%%
import xarray as xr
import numpy as np
import pandas as pd
import src.blocking.block_days as block_days
import glob
import os
import sys

#%%
import importlib
importlib.reload(block_days)
# %%
def decadal_bDays(anomaly = False):
    if anomaly:
        to_name = "bDays_anomaly.nc"
    else:
        to_name = "bDays.nc"
    print("Reading IB index...")
    IB_index = block_days.read_IB_index(anomaly = anomaly)
    print("Counting block days...")
    dec_IB = block_days.decadal_block_days(IB_index)
    print("Saving to netcdf...")
    dec_IB.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/blocking_days"+to_name)
# %%
def annual_bDays(month = 'Jun'):
    odir = f'/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_{month}/'
    files = glob.glob(odir + "*.nc")
    files.sort()
    for file in files:
        fname = os.path.basename(file)
        print(fname)
        ds = xr.open_dataset(file)
        IB_index = ds['IB index']
        ann_IB = block_days.annual_block_days(IB_index,month = month[-3:].upper())
        to_path = f'/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/bDays_{month}/{fname}'
        try:
            ann_IB.to_netcdf(to_path)
        except PermissionError:
            os.remove(to_path)
            ann_IB.to_netcdf(to_path)
            

# %%
annual_bDays(month = 'Jun')
annual_bDays(month = 'Jul')
annual_bDays(month = 'Aug')
# %%
annual_bDays(month = 'ano_Jun')
#%%
annual_bDays(month = 'ano_Jul')
annual_bDays(month = 'ano_Aug')
# %%
