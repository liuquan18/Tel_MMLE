# %%
import xarray as xr
import numpy as np
import pandas as pd

import src.composite.composite as composite
import src.reanalysis.remove_force as rmf

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
#%%
def time_convert(data):
    data['time'] = pd.to_datetime(data['time'].values)
    return data

def read_temp_data(model,var_name = 't2max'):

    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + f"{var_name}_" + month + "/"
        all_files = glob.glob(zg_path + "*.nc")
        all_files.sort()

        if model == 'ERA5' and var_name == 'ts':
            data_month = xr.open_mfdataset(all_files,combine = 'by_coords',preprocess=time_convert)
            try:
                data_month = data_month['T2M']
            except KeyError:
                data_month = data_month['var167']
        elif model == 'ERA5' and var_name == 't2max':
            data_month = xr.open_dataset(all_files[0])
            data_month = data_month['var167']
        elif model == 'ERA5' and var_name == 't2min':
            data_month = xr.open_dataset(all_files[0])
            data_month = data_month['var167']
        elif model == 'CR20':
            data_month = xr.open_dataset(all_files[0])
            data_month = data_month['air']
        elif model == 'CR20_allens':
            data_month = xr.open_mfdataset(all_files,combine = 'nested',concat_dim = 'ens')
            data_month = data_month['TMP']

        data_JJA.append(data_month)
    data = xr.concat(data_JJA, dim="time").sortby("time")
    data = data.sel(time = slice('1850','2015'))
    data_inter = rmf.detrend(data, method="quadratic_trend")
    return data_inter


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
    sel_data = composite.sel_data(period_data, extr_index)
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

def read_eof(model,group_size = 40):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    first_eof_path = odir + f"EOF_result/first_{str(group_size)}_eof_std.nc"
    last_eof_path = odir + f"EOF_result/last_{str(group_size)}_eof_std.nc"

    first_eof = xr.open_dataset(first_eof_path)
    last_eof = xr.open_dataset(last_eof_path)
    return first_eof.pc, last_eof.pc

# %%
model = "CR20_allens"
# read gph data
odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
save_path = odir + "composite/"

# %%
data = read_temp_data(model,var_name = 'ts')
First_index, Last_index = read_eof(model,group_size = 40)

#%%
n_resamples = 1000

sel_first_pos, sel_first_neg, sel_last_pos, sel_last_neg = select_data(data,First_index,Last_index)

res_first_pos, res_first_neg, res_last_pos, res_last_neg = resample_com_ind(sel_first_pos, sel_first_neg, sel_last_pos, sel_last_neg, n_resamples=n_resamples)

# %%
list_all_pros = [0] * npro
for nn in range(npro):
    list_all_pros[nn] = np.arange(n_resamples)[nn::npro]
steps = list_all_pros[rank]

#%%
for kk, step in enumerate(steps) :
    print(f"node {node}: core:{rank} kk = {kk+1}/{steps.shape[0]}")

    first_pos = sel_first_pos.isel(com=res_first_pos[:, step]).mean(dim = 'com')
    first_neg = sel_first_neg.isel(com=res_first_neg[:, step]).mean(dim = 'com')

    last_pos = sel_last_pos.isel(com=res_last_pos[:, step]).mean(dim = 'com')
    last_neg = sel_last_neg.isel(com=res_last_neg[:, step]).mean(dim = 'com')
    
    num_info = step
    print(save_path + f"first_boot/composite_first_pos_{num_info}.nc")
# %%
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
