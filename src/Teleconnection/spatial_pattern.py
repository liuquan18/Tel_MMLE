############ imports ###############
import numpy as np
import pandas as pd
import xarray as xr
from eofs.standard import Eof

import src.Teleconnection.tools as tools


def doeof(
    data: xr.DataArray,
    nmode: int = 2,
    dim: str = "com",
    standard: str = 'eof_spatial_std',
):
    """
    do eof to seasonal data along a combined dim

    **Arguments**:
        *seasondata*: The data to be decomposed, where the first dim should be the dim of 'com' or 'time'.
        *nmode*: how many modes, mode=2,means NAO and EA respectively.
        *dim*: along which dim to do the eof (calculate the co variance)
        *standard*: standard the output pc or not.
        *shuffle* random order the com (ens and time) or not.

    **Returns**:
        eof: DataArray, spatial patterns scaled (multiplied) with the temporal std of seasonal index.
             has the same spatial size as the input seasondata.[mode, lat,lon,...]
        pc: DataArray, seasonal index, if standard = True, scaled (divided) by the temporal std of itself. [mode,time]]
        exp_var: explained variance of each mode. [mode]

    """

    # make sure that there are only three dims
    data = data.squeeze()

    # make sure that the first dim is the 'com' or 'time'.
    try:
        data = data.transpose(dim, ...)
    except ValueError:
        print("no combined dimension found. use tools.stackens() first")

    # weights
    wgts = tools.sqrtcoslat(data) # same coords as data

    # EOF decompose
    solver = Eof(
        data.values, weights=wgts, center=True
    )  
    
    eof = solver.eofs(neofs=nmode)  # (mode,lat,lon,...)
    pc = solver.pcs(npcs=nmode)  # (com,mode)
    fra = solver.varianceFraction(nmode)  # (mode)

    # eof to xarray
    eofx, pcx, frax = eofs_to_xarray(data, eof, pc, fra)

    # deweight
    eofx = eofx / wgts[0]

    # standardize, here the loading gives to the pc, to make the index from different spatil pattern comparable.
    if standard == "eof_spatial_std":
        eofx, pcx = standard_by_eof_spatial_std(eofx, pcx)
    elif standard == "pc_temporal_std":
        eofx, pcx = standard_by_pc_temporal_std(eofx, pcx)

    # fix the sign, so that the North center of action is always low.
    eofx, pcx = fix_sign(eofx, pcx)

    # make sure at the loc where the data==np.nan, the eof==np.nan as well.
    map_data = data[0]  # just one map
    eofx = eofx.where(np.logical_not(map_data.isnull()), map_data)

    # unstack the dim 'ens' and 'time' or 'win'
    pcx = pcx.unstack()

    # dorp vars
    try:
        eofx = eofx.drop_vars(("ens", "time", dim))
    except ValueError:
        pass

    # to dataset
    eof_result = xr.Dataset({"eof": eofx, "pc": pcx, "fra": frax})

    return eof_result

def standard_by_eof_spatial_std(eofx, pcx):
    std_eof = eofx.std(dim=("lat", "lon"))
    eofx = eofx / std_eof
    std_eof = std_eof.squeeze()
    pcx = pcx * std_eof
    return eofx,pcx

def standard_by_pc_temporal_std(eofx, pcx):
    std_pc = pcx.std(dim="time")
    pcx = pcx / std_pc
    std_pc = std_pc.squeeze()
    eofx = eofx * std_pc
    return eofx,pcx

def eofs_to_xarray(data, eof, pc, fra):
    reduce_dim = data.dims[0] # 'com' or 'time'
    eof_cnt = data[0]
    time_tag = data.unstack().time.values[0] # the time tag of the first map
    eof_cnt = eof_cnt.drop_vars(reduce_dim) # drop the dim 'com' or 'time'
    eof_cnt = eof_cnt.expand_dims(dim = {'mode':['NAO','EA'],'decade':[time_tag]},axis = [0,-1]) # add one more dimension to eof_cnt with shape 1

    eof = eof[..., np.newaxis] # add one new dimension for the info of decade (time of the beginning of the decade)
    eofx = eof_cnt.copy(data=eof)

    # pc to xarray
    pcx = xr.DataArray(
        pc, dims=[reduce_dim, "mode"], coords={reduce_dim: data[reduce_dim], "mode": ["NAO", "EA"]}
    )
    frax = xr.DataArray(fra, dims=["mode"], coords={"mode": ["NAO", "EA"]})
    return eofx,pcx,frax


def fix_sign(eof,pc):
    """
    function to calculate the coefficient for eof, so that the sign is consistent.
    for NAO, the positive NAO with a low in the North and high at the south.
    for EA, the positive EA with a low in the center.
    **Arguments**:
        eof: the eof (spatial pattern) to be changed. much include NAO and EA the same time.
    **Returns**:
        coefficient of NAO and EA in xarray.
    """

    # sortby lat since some dataset the lat goes from higher to lower.
    eof = eof.sortby("lat")

    # NAO
    coef_NAO = (
        eof.sel(lat=slice(60, 90), lon=slice(-70, -10), mode="NAO").mean(
            dim=["lat", "lon"]
        )
        < 0
    )
    coef_NAO = 2 * coef_NAO - 1  # to make 1 to 1 , 0 to -1

    # EA
    coef_EA = (
        eof.sel(lat=slice(45, 65), lon=slice(-40, 40), mode="EA").mean(
            dim=["lat", "lon"]
        )
        < 0
    )
    coef_EA = 2 * coef_EA - 1

    # change sign
    coef = xr.concat([coef_NAO, coef_EA], dim="mode")
    eofx = eofx * coef
    coef = coef.squeeze()
    pcx = pcx * coef
    return eofx, pcx
