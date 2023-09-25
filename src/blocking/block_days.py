#%%
import xarray as xr
import numpy as np
import pandas as pd


# %%
def read_IB_index(model = 'MPI_GE_onepct_30_daily',anomaly = False):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/"
    if anomaly:
        months = ['ano_Jun','ano_Jul','ano_Aug']
    else:
        months = ['Jun','Jul','Aug']
    blocks = []
    for month in months:
        fs = odir + f"block_{month}/*.nc"
        ds = xr.open_mfdataset(fs, combine='nested', concat_dim='ens')
        block = ds['IB index']
        blocks.append(block)
    blocks = xr.concat(blocks, dim = 'time')
    blocks = blocks.sortby('time')
    return blocks
# %%
def decadal_block_days(IB_index):
    dec_IB = IB_index.resample(time = '10AS-JUN').sum(dim = ('time','ens'))
    dec_IB = dec_IB[:-1] # exclude the last decade where there is no full decade
    return dec_IB
def annual_block_days(IB_index):
    ann_IB = IB_index.resample(time = 'AS-JUN').sum(dim = ('time'))
    return ann_IB

