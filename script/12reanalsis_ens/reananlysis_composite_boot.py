# %%
import xarray as xr
import numpy as np
import pandas as pd

import src.composite.composite as composite

import mpi4py as MPI
import glob
import os
import sys

# %%
import importlib

importlib.reload(composite)


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
node = int(sys.argv[1])
f1 = int(sys.argv[2])
f2 = int(sys.argv[3])

# %%
def linear_detrend(data):
    ens_mean = data.mean(dim="ens")
    ens_mean_c = ens_mean.copy()
    ens_mean_c["time"] = np.arange(ens_mean_c.time.size)
    linear_coef = ens_mean_c.polyfit(dim="time", deg=1)
    linear_fitted = xr.polyval(ens_mean_c.time, linear_coef.polyfit_coefficients)
    linear_fitted["time"] = ens_mean.time
    linear_detrend = data - linear_fitted
    return linear_detrend


def read_ts_data(model, external_forcing="linear_trend"):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + "ts_" + month + "/"
        data_month = xr.open_mfdataset(
            zg_path + "*.nc", combine="nested", concat_dim="ens"
        )["ts"]
        data_detrend = linear_detrend(data_month)
        data_JJA.append(data_detrend)
    data = xr.concat(data_JJA, dim="time").sortby("time")

    data = data.chunk(chunks={"time": 100})
    return data


# %%
def _select_data(data, extr_type, index):
    period_data = data.sel(time=index.time, method="nearest")
    period_data = period_data.sortby("time")
    period_data["time"] = index.time


    # get the coordinates of the extremes
    extr_index = composite.extreme(index, extreme_type=extr_type, threshold=1.5)
    extr_index.attrs["extreme_type"] = extr_type
    extr_index = extr_index.sel(mode="NAO")

    extr_index = extr_index.sortby("time")
    sel_data = composite._sel_data(period_data, extr_index)
    return sel_data

def select_data(data,First_index,Last_index):
    sel_first_pos = _select_data(data, "pos", First_index)
    sel_first_neg = _select_data(data, "neg", First_index)

    sel_last_pos = _select_data(data, "pos", Last_index)
    sel_last_neg = _select_data(data, "neg", Last_index)

    return sel_first_pos, sel_first_neg, sel_last_pos, sel_last_neg
#%%
def _resample_com_ind(sel_data, n_resamples=1000):
    n_samples = sel_data.com.size
    rng = np.random.default_rng(seed=12345)
    sampled_index = rng.choice(n_samples, size=(n_samples, n_resamples), replace=True)
    return sampled_index

def resample_com_ind(sel_first_pos, sel_first_neg, sel_last_pos, sel_last_neg,n_resamples=1000):
    res_first_pos = _resample_com_ind(sel_first_pos, n_resamples=n_resamples)
    res_first_neg = _resample_com_ind(sel_first_neg, n_resamples=n_resamples)

    res_last_pos = _resample_com_ind(sel_last_pos, n_resamples=n_resamples)
    res_last_neg = _resample_com_ind(sel_last_neg, n_resamples=n_resamples)
    return res_first_pos, res_first_neg, res_last_pos, res_last_neg
        
# %%
model = "ERA5_allens"
# read gph data
odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
save_path = odir + "composite/"

# %%
data = read_ts_data("ERA5_allens")
First_index = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/EOF_result/first_40_eof_std.nc"
).pc
Last_index = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/EOF_result/last_40_eof_std.nc"
).pc

#%%
n_resamples = 50

sel_first_pos, sel_first_neg, sel_last_pos, sel_last_neg = select_data(data,First_index,Last_index)

res_first_pos, res_first_neg, res_last_pos, res_last_neg = resample_com_ind(sel_first_pos, sel_first_neg, sel_last_pos, sel_last_neg, n_resamples=n_resamples)

# %%
list_all_pros = [0] * npro
for nn in range(npro):
    list_all_pros[nn] = np.arange(n_resamples)[nn::npro]
steps = list_all_pros[rank]

#%%
for kk, step in enumerate(steps) :
    print(f"node {node}: kk = {kk+1}/{steps.shape[0]}")

    first_pos = sel_first_pos.isel(com=res_first_pos[:, step]).mean(dim = 'com')
    first_neg = sel_first_neg.isel(com=res_first_neg[:, step]).mean(dim = 'com')

    last_pos = sel_last_pos.isel(com=res_last_pos[:, step]).mean(dim = 'com')
    last_neg = sel_last_neg.isel(com=res_last_neg[:, step]).mean(dim = 'com')
    
    num_info = node * n_resamples + step

    print("saving files ...")
    first_pos.to_netcdf(
        save_path + f"first_boot/composite_first_pos_{num_info}.nc"
    )

    first_neg.to_netcdf(
        save_path + f"first_boot/composite_first_neg_{num_info}.nc"
    )

    last_pos.to_netcdf(
        save_path + f"last_boot/composite_last_pos_{num_info}.nc"
    )

    last_neg.to_netcdf(
        save_path + f"last_boot/composite_last_neg_{num_info}.nc"
    )
# %%
