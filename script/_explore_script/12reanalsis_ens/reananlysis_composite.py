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
client, cluster = slurm_cluster.init_dask_slurm_cluster(scale = 2, processes = 1, memory = "512GB", walltime = "08:00:00")


#%%
def time_convert(data):
    data['time'] = pd.to_datetime(data['time'].values)
    return data

def read_temp_data(model,var_name = 't2max', var_code = None):
    if var_code is None:
        var_code = var_name

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
            try:
                data_month = data_month['TMP']
            except KeyError:
                data_month = data_month[var_code]

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
def composite_reana(model,var_name = 'ts', group_size = 40,reduction = 'mean'):
    var_data = read_temp_data(model,var_name = var_name)
    first_index, last_index = read_eof(model,group_size=group_size)
    first_composite = composite.Tel_field_composite(
        first_index,
        var_data,
        threshold=1.5,
        reduction=reduction,
        bootstrap=False,
    )
    last_composite = composite.Tel_field_composite(
        last_index,
        var_data,
        threshold=1.5,
        reduction=reduction,
        bootstrap=False,
    )
    diff = last_composite - first_composite
    return first_composite, last_composite, diff


#%%
def composite_oneclick(model, var_name,group_size,reduction = 'mean',save = False):
    ERA_first, ERA_last, ERA_diff = composite_reana(model,var_name = var_name, group_size=group_size,reduction = reduction)
    ERA_first.name = 'ts'
    ERA_last.name = 'ts'
    ERA_diff.name = 'ts'
    # composite_plot.composite_plot(ERA_first, ERA_last, 'NAO', levels=np.arange(-1.5, 1.6, 0.3))

    if save:
        ERA_first.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/composite/first_composite_{reduction}_{var_name}_{group_size}.nc")
        ERA_last.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/composite/last_composite_{reduction}_{var_name}_{group_size}.nc")
    return ERA_first, ERA_last, ERA_diff
#%%

def composite_together(model,var_name,group_size,reduction = 'mean',alpha = 0.05):
    """
    first,last,diff,first_boot,last_boot,diff_sig
    """
    odir = f'/work/mh0033/m300883/Tel_MMLE/data/{model}/composite/'
    first_composite = xr.open_dataset(odir + f'first_composite_{reduction}_{var_name}_{group_size}.nc')
    last_composite = xr.open_dataset(odir + f'last_composite_{reduction}_{var_name}_{group_size}.nc')
    diff = last_composite - first_composite

    first_composite_boot_pos = xr.open_mfdataset(
        odir + 'first_boot/composite_first_pos*.nc',
        combine = 'nested',
        concat_dim = 'bootstrap',
    ).chunk(dict(bootstrap=-1))

    first_composite_boot_neg = xr.open_mfdataset(
        odir + 'first_boot/composite_first_neg*.nc',
        combine = 'nested',
        concat_dim = 'bootstrap',
    ).chunk(dict(bootstrap=-1))
    first_composite_boot = xr.concat([first_composite_boot_pos['__xarray_dataarray_variable__'],
                                      first_composite_boot_neg['__xarray_dataarray_variable__']],
                                      pd.Index(['pos','neg'],name = 'extr_type'))


    last_composite_boot_pos = xr.open_mfdataset(
        odir + 'last_boot/composite_last_pos*.nc',
        combine = 'nested',
        concat_dim = 'bootstrap',
    ).chunk(dict(bootstrap=-1))
    last_composite_boot_neg = xr.open_mfdataset(
        odir + 'last_boot/composite_last_neg*.nc',
        combine = 'nested',
        concat_dim = 'bootstrap',
    ).chunk(dict(bootstrap=-1))
    last_composite_boot = xr.concat([last_composite_boot_pos['__xarray_dataarray_variable__'],
                                     last_composite_boot_neg['__xarray_dataarray_variable__']],
                                     pd.Index(['pos','neg'],name = 'extr_type'))

    diff_boot = last_composite_boot - first_composite_boot

    # check if the difference is significant
    # get the 95% confidence interval
    low_bnd = alpha/2.0
    high_bnd = 1-alpha/2.0
    ci = diff_boot.quantile([low_bnd, high_bnd], dim="bootstrap")
    # check if 0 is in the interval, return boolean
    diff_sig = xr.where((ci.sel(quantile=low_bnd) > 0) | (ci.sel(quantile=high_bnd) < 0), 1, 0)
    diff_sig.name = var_name
    # combine the first, last, diff and diff_sig together
    period = xr.IndexVariable('period', ['first', 'last', 'diff', 'diff_sig'])
    composite = xr.concat([first_composite['ts'], 
                           last_composite['ts'], 
                           diff['ts'], 
                           diff_sig], 
                           dim=period,
                           coords='minimal')

    # save
    composite.to_netcdf(odir + f'composite_{reduction}_{var_name}_{group_size}.nc')
# %%
composite_reana('CR20_allens', 'psl')