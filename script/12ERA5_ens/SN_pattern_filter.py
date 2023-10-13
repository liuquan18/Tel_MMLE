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
from eofs.standard import Eof

#%%
import importlib
importlib.reload(index_generate)

#%%
client, cluster = slurm_cluster.init_dask_slurm_cluster(scale = 2, processes = 1, memory = "512GB", walltime = "04:30:00")


#%%
model = 'ERA5_allens'
var_name = "zg"
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
    zg_path = odir + f"{var_name}_" + month + "/"
    data_month = index_generate.read_data(zg_path, plev=plev, remove_ens_mean=False,var_name=var_name)
    data_JJA.append(data_month)
data = xr.concat(data_JJA, dim="time").sortby("time")
#%%
data = data.sel(time = slice("1940","2022"))

#%%
# data = data.coarsen(lat=2, lon=2, boundary="trim").mean()

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
# %%
# keep unscaled copies of these variables
Xt_ensmean = data_anomaly.mean("ens")
Xt_flat = data_anomaly.stack(index=["time", "ens"]).stack(shape=["lat", "lon"])
Xt_ensmean_flat = Xt_ensmean.stack(shape=["lat", "lon"])

index = X_flat.index
n = len(index)
# %%
    # Large Ensemble EOFs
# def ens_eof(X_flat, n):
#     u,s,v = np.linalg.svd(np.transpose(X_flat.values)/np.sqrt(n-1))
#     return s, v

def ens_eof(X_flat,neof = 200):
    solver = Eof(X_flat.values,center=False)
    pcvec = solver.eofs(neofs=neof)
    pcvec = pcvec.T # transpose to match the np.linalg.eig output
    evl = solver.eigenvalues(neigs=neof)
    s = np.sqrt(evl)
    return s, pcvec

#%%
# Large Ensemble Forced Patterns
def SN_pattern(s,pcvec, neof = 200):
    S = np.matmul(pcvec[:, 0:neof], np.diag(1 / s[0:neof]))
    Sadj = np.matmul(np.diag(s[0:neof]), pcvec[:, 0:neof].T)

    ensmeanPCs = np.matmul(X_ensmean_flat.values, S)  # ensemble-mean principal components

    gamma = np.cov(ensmeanPCs.T)  # covariance matrix of ensemble-mean principal components

    lat = X.lat
    lon = X.lon
    scale = wgts.isel(ens = 0, time = 0,lon = 0)

    u2, signal_frac, v2 = np.linalg.svd(gamma)
    SNP = np.matmul(v2, Sadj)
    SNPs_reshaped = SNP.reshape(neof, len(lat), len(lon)) / scale.values[None, :, None]

    weights = np.matmul(S, v2.T)
    weights = weights.reshape(len(lat), len(lon), neof) * scale.values[:, None, None]
    weights = weights.reshape(len(lat) * len(lon), neof)

    tk = np.matmul(Xt_flat.values, weights)  # compute timeseries from full data matrix

    tk_emean = np.matmul(
        Xt_ensmean_flat.values, weights
    )  # compute ensemble-mean timeseries from ensemble-mean data

    return SNPs_reshaped,tk_emean,tk,signal_frac
# %%
s, pcvec = ens_eof(X_flat, neof = 200)
# # %%
# import pickle
# outputdir = "/work/mh0033/m300883/Tel_MMLE/script/12ERA5_ens/ERA5_tickle/" 
# pickle.dump(
#     [s, v],
#     open(outputdir + "pyeofs_ensmean.pkl", "wb"),
#     protocol=4,
# )
# %%
neof = 200  # number of EOFs retained in S/N maximizing pattern analysis
SNPs_reshaped,tk_emean,tk,signal_frac = SN_pattern(s, pcvec, neof = neof)
# %%
M = 50
lat = X.lat
lon = X.lon
T = X.time
scale = wgts.isel(ens = 0, time = 0,lon = 0)
# %%
X_forced = np.matmul(
    tk_emean[:, 0:M], SNPs_reshaped[0:M, :, :].reshape(M, len(lat) * len(lon))
)
#%%
X_forced = X_forced.reshape(len(T), len(lat), len(lon))
# %%
X_forced = X.isel(ens = 0).copy(data = X_forced)

#%%
X_forced.to_netcdf(save_path + f"X_forced_{var_name}_200_50.nc")
data_anomaly.to_netcdf(save_path + f"X_anomaly_{var_name}.nc")

# %%
internal = data_anomaly - X_forced
# %%
def decompose(xarr):
    field = xarr.sortby('time')
    field = field.stack(com = ('ens','time'))
    eof_result = ssp.doeof(field,standard = 'eof_spatial_std')
    return eof_result

# %%
first_40 = internal.sel(time=slice('1940','1980'))

last_40 = internal.sel(time=slice('1982','2022'))

first_40_eof = decompose(first_40)
last_40_eof = decompose(last_40)

# %%
def standard_40(first_40_eof,last_40_eof):
    mean = first_40_eof['pc'].mean(dim=('time','ens'))
    std = first_40_eof['pc'].std(dim=('time','ens'))

    first_40_eof['pc'] = (first_40_eof['pc'] - mean) / std
    last_40_eof['pc'] = (last_40_eof['pc'] - mean) / std
    return first_40_eof, last_40_eof

#%%

first_40_eof_std, last_40_eof_std = standard_40(first_40_eof,last_40_eof)
#%%

fig2 = pplt.figure(figsize=(180 / 25.4, 90 / 25.4),sharex=False,sharey=False)
fig2.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
    
)
gs = pplt.GridSpec(
    ncols=2,
    nrows=1,
)
ax1 = fig2.add_subplot(gs[0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax2 = fig2.add_subplot(gs[1])


ax1, fmap, lmap = stat_overview.spatial_pattern_plot(
    ax1,
    first_40_eof_std.eof.sel(mode = 'NAO',decade = '1940').squeeze(),
    first_40_eof_std.fra.sel(mode = 'NAO',decade = '1940').squeeze(),
    last_40_eof_std.eof.sel(mode = 'NAO',decade = '1982').squeeze(),
    last_40_eof_std.fra.sel(mode = 'NAO',decade = '1982').squeeze(),
    levels = np.arange(-2,2.1,0.4),
)
ax2,hist = stat_overview.index_distribution_plot(
    ax2,
    first_40_eof_std.pc.sel(mode = 'NAO'),
    last_40_eof_std.pc.sel(mode = 'NAO'),
)
# %%
