#  %%
import xarray as xr
import numpy as np
import pandas as pd
import random
import os

import src.warming_stage.warming_stage as warming_stage
import src.Teleconnection.spatial_pattern as ssp
import src.MMLE_TEL.index_generator as index_generate

import glob

import proplot as pplt
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import matplotlib.pyplot as plt
import src.Teleconnection.tools as tools
import src.compute.slurm_cluster as slurm_cluster
import pickle

#%%
import importlib
importlib.reload(index_generate)

#%%
client, cluster = slurm_cluster.init_dask_slurm_cluster(scale = 2, processes = 1, memory = "512GB", walltime = "04:30:00")

#%%
model = 'ERA5_allens'
fixedPattern="decade"
plev=50000
standard="first"
#%%
# read gph data
odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
ts_mean_path = odir + "ts_processed/ens_fld_year_mean.nc"
save_path = odir + "EOF_result/"

#%%
data_JJA = []
for month in ["Jun", "Jul", "Aug"]:
    print(f"reading the gph data of {month} ...")
    zg_path = odir + "zg_" + month + "/"
    data_month = index_generate.read_data(zg_path, plev=plev, remove_ens_mean=False)
    data_JJA.append(data_month)
data = xr.concat(data_JJA, dim="time").sortby("time")
#%%
data = data.fillna(0)
#%%
data_anomaly = data.groupby('time.month') - data.groupby('time.month').mean('time')
#%%
# pre_process the data
# weights
wgts = tools.sqrtcoslat(data) # same coords as data

# %%
X = data_anomaly * wgts
X_ensmean = X.mean(dim='ens')
X_flat = X.stack(index=('ens','time')).stack(shape=('lat','lon'))
X_ensmean_flat = X_ensmean.stack(shape=('lat','lon'))

# #%%
# # rechunk
# X_flat = X_flat.chunk({'index': 2520, 'shape': 1000})
# X_ensmean_flat = X_ensmean_flat.chunk({'shape': 1000})


# %%
# keep unscaled copies of these variables
Xt_ensmean = data_anomaly.mean("ens")
Xt_flat = data_anomaly.stack(index=["time", "ens"]).stack(shape=["lat", "lon"])
Xt_ensmean_flat = Xt_ensmean.stack(shape=["lat", "lon"])

# # rechunk
# Xt_flat = Xt_flat.chunk({'index':  2520, 'shape': 1000})
# Xt_ensmean_flat = Xt_ensmean_flat.chunk({'shape': 1000})

#%%

index = X_flat.index
n = len(index)
# %%
# Large Ensemble EOFs
def ens_eof(X_flat, n):
    Cov = np.matmul(X_flat.values.T, X_flat.values) / (n - 1)
    evl, pcvec = np.linalg.eig(Cov)
    s = np.sqrt(evl)
    return s, pcvec
# %%
s, pcvec = ens_eof(X_flat, n)
# %%
import pickle
outputdir = "/work/mh0033/m300883/Tel_MMLE/script/12ERA5_ens/ERA5_tickle/" 
pickle.dump(
    [s, pcvec],
    open(outputdir + "EOF_ensmean.pkl", "wb"),
    protocol=4,
)

# %%
