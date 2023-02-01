#%% 
# import xarray, numpy, pandas and matplotlib
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import src.extreme.period_pattern_extreme as extreme


#%%
dir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/ind_first_pc.nc"
ds = xr.open_dataset(dir)

# %%
# change ds time to datetime
ds['time'] = ds.indexes['time'].to_datetimeindex()

# normalize the data
ds = (ds - ds.mean(dim='time')) / ds.std(dim='time')

pc = ds.pc
# %%
# select the first 10 years of pc
pc_first = pc.sel(time = slice('1850-01-01', '1860-01-01'))
# %%
count = extreme.period_extreme_count(pc_first,pc)
# %%
