# %%
import xarray as xr
import numpy as np
import pandas as pd
import random
import os

import src.warming_stage.warming_stage as warming_stage
import src.Teleconnection.spatial_pattern as ssp
import src.MMLE_TEL.index_generator as index_generate
import src.composite.composite as composite
import glob

import proplot as pplt
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import matplotlib.pyplot as plt
#%%
import src.compute.slurm_cluster as slurm_cluster
import pickle

#%%
import importlib
importlib.reload(composite)

#%%
client, cluster = slurm_cluster.init_dask_slurm_cluster(scale = 2, processes = 1, memory = "512GB", walltime = "08:00:00")


#%%
def linear_detrend(data):
    ens_mean = data.mean(dim = 'ens')
    ens_mean_c = ens_mean.copy()
    ens_mean_c['time'] = np.arange(ens_mean_c.time.size)
    linear_coef = ens_mean_c.polyfit(dim='time', deg=1)
    linear_fitted = xr.polyval(ens_mean_c.time, linear_coef.polyfit_coefficients)
    linear_fitted['time'] = ens_mean.time
    linear_detrend = data - linear_fitted
    return linear_detrend

#%%
model = 'ERA5_allens'

#%%
# read gph data
odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
save_path = odir + "composite/"

data_JJA = []
for month in ["Jun", "Jul", "Aug"]:
    print(f"reading the gph data of {month} ...")
    zg_path = odir + "ts_" + month + "/"
    data_month = xr.open_mfdataset(zg_path + "*.nc",combine = 'nested',concat_dim = 'ens')['ts']
    data_detrend = linear_detrend(data_month)
    data_JJA.append(data_detrend)
data = xr.concat(data_JJA, dim="time").sortby("time")
# rechunk
#%%
# data = data.chunk(chunks = {'time': 100})

# %%
First_index = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/EOF_result/first_40_eof_std.nc").pc
Last_index = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/EOF_result/last_40_eof_std.nc").pc

# %%

# first 40 years
first_composite = composite.Tel_field_composite(
    First_index,
    data,
    threshold=1.5,
    reduction='mean',
    bootstrap=False,
)

# %%
last_composite = composite.Tel_field_composite(
    Last_index,
    data,
    threshold=1.5,
    reduction='mean',
    bootstrap=False,
)
#%%
first_composite.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/composite/first_composite.nc")
last_composite.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/composite/last_composite.nc")


# %%
# combine the first and last composite
period = xr.IndexVariable('period', ['first', 'last'])
com_ts = xr.concat([first_composite, last_composite], dim=period)
diff = last_composite - first_composite

# %%
print(" doing bootstrap resampling...")
# first 10 years with bootstrap

# first 40 years
first_composite_boot = composite.Tel_field_composite(
    First_index,
    data,
    threshold=1.5,
    reduction='mean',
    bootstrap=True,
)

# %%
last_composite_boot = composite.Tel_field_composite(
    Last_index,
    data,
    threshold=1.5,
    reduction='mean',
    bootstrap=True,
)
#%%
first_composite_boot.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/composite/first_composite_boot.nc")
#%%
last_composite_boot.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/composite/last_composite_boot.nc")

#%%
# difference between first and last 10 years
diff_boot = last_composite_boot - first_composite_boot
#%%
# check if the difference is significant
# get the 95% confidence interval
alpha = 0.05
low_bnd = alpha/2.0
high_bnd = 1-alpha/2.0

#%%
diff_boot.compute()

#%%
ci = diff_boot.quantile([low_bnd, high_bnd], dim="bootstrap")
# check if 0 is in the interval, return boolean
diff_sig = xr.where((ci.sel(quantile=low_bnd) > 0) | (ci.sel(quantile=high_bnd) < 0), 1, 0)
# combine the first, last, diff and diff_sig together
period = xr.IndexVariable('period', ['first', 'last', 'diff', 'diff_sig'])
Com_ts = xr.concat([first_composite, last_composite, diff, diff_sig], dim=period)
# %%
