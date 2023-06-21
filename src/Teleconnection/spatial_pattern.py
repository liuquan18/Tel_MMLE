############ imports ###############
import numpy as np
import pandas as pd
import xarray as xr
from eofs.standard import Eof
from tqdm.notebook import tqdm, trange

import src.Teleconnection.tools as tools


def doeof(
    data: xr.DataArray,
    nmode: int = 2,
    dim: str = "com",
    standard: bool = True,
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
    wgts = tools.sqrtcoslat(data)

    # EOF decompose
    solver = Eof(
        data.values, weights=wgts, center=True
    )  # if it's com dim, is is right to remove the mean
    # along the com dim?
    eof = solver.eofs(neofs=nmode)  # (mode,lat,lon,...)
    pc = solver.pcs(npcs=nmode)  # (com,mode)
    fra = solver.varianceFraction(nmode)  # (mode)

    # eof to xarray
    eof_cnt = data.unstack()
    eof_cnt = eof_cnt.isel(ens=[0, 1], time=[0])
    eof_cnt = eof_cnt.rename({"ens": "mode", "time": "decade"})
    eof_cnt = eof_cnt.transpose("mode", ...)
    eof_cnt["mode"] = ["NAO", "EA"]

    eof = eof[..., np.newaxis]
    eofx = eof_cnt.copy(data=eof)

    # pc to xarray
    pcx = xr.DataArray(
        pc, dims=[dim, "mode"], coords={dim: data[dim], "mode": ["NAO", "EA"]}
    )
    frax = xr.DataArray(fra, dims=["mode"], coords={"mode": ["NAO", "EA"]})

    # deweight
    eofx = eofx / wgts.isel(com=0)

    # standardize, here the loading gives to the pc, to make the index from different spatil pattern comparable.
    std_eof = eofx.std(dim=("lat", "lon"))
    eofx = eofx / std_eof
    std_eof = std_eof.squeeze()
    pcx = pcx * std_eof

    # change sign
    coef = sign_coef(eofx)
    eofx = eofx * coef
    coef = coef.squeeze()
    pcx = pcx * coef

    # make sure at the loc where the data==np.nan, the eof==np.nan as well.
    map_data = data[0]  # just one map
    eofx = eofx.where(np.logical_not(map_data.isnull()), map_data)

    # unstack the dim 'ens' and 'time' or 'win'
    pcx = pcx.unstack()

    # dorp vars
    eofx = eofx.drop_vars(("ens", "time", "com"))

    # to dataset
    eof_result = xr.Dataset({"eof": eofx, "pc": pcx, "fra": frax})

    return eof_result


def sign_coef(eof):
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

    return xr.concat([coef_NAO, coef_EA], dim="mode")
