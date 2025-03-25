#%%
import xarray as xr
import numpy as np
import pandas as pd
import glob
import warnings
from src.MMLE_TEL.index_stats import decadal_extrc
#%%

def read_data(
    zg_path,
    plev=None,
    remove_ens_mean=True,
    var_name = "zg",
):
    """
    read data quickly
    """
    gph_dir = zg_path
    # read MPI_onepct data
    # fix the order of ensemble members
    print("reading the gph data of all ensemble members...")
    all_ens_lists = sorted(
        glob.glob(gph_dir + "*.nc")
    )  # to make sure that the order of ensemble members is fixed
    zg_data = xr.open_mfdataset(
        all_ens_lists, combine="nested", concat_dim="ens", join="override",
    )  # consider chunks={}, # kz the file size is small (< 3G). 

    zg_data["ens"] = np.arange(zg_data.ens.size)
    try:
        zg_data = zg_data.var156
    except AttributeError:
        try:
            zg_data = zg_data.zg
    
        except AttributeError:
            zg_data = zg_data[var_name]

    # time to datetime
    try:
        zg_data["time"] = zg_data.indexes["time"].to_datetimeindex()
    except AttributeError:
        zg_data["time"] = pd.to_datetime(zg_data.time)

    # demean
    if remove_ens_mean:
        print(" demean the ensemble mean...")
        zg_ens_mean = zg_data.mean(dim="ens")
        zg_demean = zg_data - zg_ens_mean
    else:
        zg_demean = zg_data

    # select one altitude
    try:
        if plev is not None:
            print(" select the specific plev...")
            zg_plev = zg_demean.sel(plev=plev)
        else:
            # select the 1000hPa - 200hPa
            print(" select the 1000hPa - 200hPa...")
            zg_plev = zg_demean.sel(plev=slice(100000, 20000))
            if zg_plev.plev.size == 0:
                zg_plev = zg_demean.sel(plev=slice(20000, 100000))
    except KeyError:
        zg_plev = zg_demean # for the data only with one plev
    return zg_plev

# %%
def box_diff(zg):
    south_box = zg.sel(lat = slice(55, 45), lon = slice(-25, 5)).mean(dim = ('lat', 'lon'))
    north_box = zg.sel(lat = slice(70, 60), lon = slice(-52, -22)).mean(dim = ('lat', 'lon'))
    diff = south_box - north_box
    return diff

# %%
# read gph data

odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/"
data_JJA = []
for month in ["Jun", "Jul", "Aug"]:
    print(f"reading the gph data of {month} ...")
    zg_path = odir + "zg_" + month + "/"
    data_JJA.append(read_data(zg_path, plev=50000))
data = xr.concat(data_JJA, dim="time").sortby("time")
# %%
data = data.sortby('time')
# %%
NAO_index = box_diff(data)
# %%
NAO_first = NAO_index.sel(time = slice('1850','1859'))
mean = NAO_first.mean()
std = NAO_first.std()

NAO_last = NAO_index.sel(time = slice('2090','2099'))
# %%
NAO_index_ano = (NAO_index - mean)/std
# %%
NAO_index_ano.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/box_diff/NAO_index_ano.nc")
# %%
extrc = decadal_extrc(NAO_index_ano)
# %%
extrc.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/box_diff/NAO_index_ano_extre_counts.nc")
# %%
