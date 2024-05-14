# %%
import xarray as xr
import numpy as np
import pandas as pd

import src.composite.composite as composite
import glob

import src.reanalysis.remove_force as remove_force
import src.plots.composite_plot as composite_plot

#%%
import src.compute.slurm_cluster as slurm_cluster
import pickle

#%%
import importlib
importlib.reload(composite)
importlib.reload(composite_plot)

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
    data_inter = remove_force.detrend(data, method="quadratic_trend")
    return data_inter

#%%
def read_eof(model,group_size = 40):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    first_eof_path = odir + f"EOF_result/first_{str(group_size)}_eof_std.nc"
    last_eof_path = odir + f"EOF_result/last_{str(group_size)}_eof_std.nc"

    first_eof = xr.open_dataset(first_eof_path)
    last_eof = xr.open_dataset(last_eof_path)
    return first_eof.pc, last_eof.pc

#%%
def composite_reana(model,var_name = 'ts', group_size = 40,reduction = 'mean',**kwargs):
    var_data = read_temp_data(model,var_name = var_name)
    var_data.load()
    first_index, last_index = read_eof(model,group_size=group_size)
    print(f" compositing {reduction} of the {var_name} data...")
    composite_mean = composite.first_last_extreme_composite(
        first_index,
        last_index,
        var_data,
        threshold=1.5,
        reduction=reduction,
        return_diff=True,
        **kwargs,
    )
    save_path = f'/work/mh0033/m300883/Tel_MMLE/data/{model}/composite/composite_{reduction}_{var_name}_{group_size}_withboot.nc'
    composite_mean.to_netcdf(save_path)

# %%
if __name__ == "__main__":
    model = 'CR20_allens'
    var_name = 'ts'
    group_size = 40
    reduction = 'mean'
    composite_reana(model,var_name = var_name, group_size = group_size,reduction = reduction)


# %%
