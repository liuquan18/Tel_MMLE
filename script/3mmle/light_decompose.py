#%%
# imports
import xarray as xr
import numpy as np
import pandas as pd

import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools

#%%
gph_dir = "/work/mh0033/m300883/Tel_MMLE/data/CanESM2/zg_processed/"
period = "first"
if period == "first":
    period = slice("1951-01-01", "1960-12-31")

#%%
def read_data(gph_dir):
    print("reading data...")
    zg_data = xr.open_mfdataset(
        gph_dir + "*.nc", combine="nested", concat_dim="ens", join="override"
    )
    try:
        zg_data = zg_data.var156
    except AttributeError:
        zg_data = zg_data.zg

    # time to datetime
    try:
        zg_data["time"] = zg_data.indexes["time"].to_datetimeindex()
    except AttributeError:
        zg_data["time"] = pd.to_datetime(zg_data.time)

    # demean
    print(" demean the ensemble mean...")
    zg_ens_mean = zg_data.mean(dim="ens")
    zg_demean = zg_data - zg_ens_mean

    # select trop
    print(" select troposphere...")
    zg_trop = zg_demean.sel(plev=slice(100000, 20000))
    if zg_trop.plev.size == 0:
        zg_trop = zg_demean.sel(plev=slice(20000, 100000))

    # standardize seperately with the temporal mean and std
    print(" standardize each altitudes seperately...")
    zg_trop = (zg_trop - zg_trop.mean(dim="time")) / zg_trop.std(dim="time")
    return zg_trop


# %%

# the spatial pattern
def spatial_pattern(zg_data, period):

    # select period
    data = zg_data.sel(time=period)

    # combine ens and time
    data = data.stack(com=("ens", "time"))

    # random order the data along the combined dimension
    data = data.sel(com=np.random.permutation(data.com))

    eof, _, fra = ssp.doeof(data, nmode=2, dim="com", standard=True)
    return eof, fra


def pc(zg_data, eofs):
    data = zg_data.stack(com=("ens", "time"))
    pc = ssp.project_field(data, eofs, dim="com")
    return pc


#%%

# spatial pattern
def decompose(zg_data, period):
    """
    decompose the 3d data into EOFs, FRAs and PCs
    parameters:
    zg_data: xarray.DataArray
        3d data
    period: slice of the time where the spatial patterns are calculated
    return:
    eofs: xarray.DataArray
        EOFs
    fras: xarray.DataArray
        FRAs
    pcs: xarray.DataArray
        PCs
    """
    print("decompose the spatial pattern...")
    eofs = zg_data.groupby("plev", squeeze=True).apply(
        lambda x, period=period: spatial_pattern(x, period)[0]
    )
    fras = zg_data.groupby("plev", squeeze=True).apply(
        lambda x, period=period: spatial_pattern(x, period)[1]
    )
    print("decompose the temporal indexes...")
    pcs = pc(zg_data, eofs)
    return eofs, fras, pcs


# %%
zg_data = read_data(gph_dir)

#%%
eofs, fras, pcs = decompose(zg_data, period)
# %%
# weight th
