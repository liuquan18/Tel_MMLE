import numpy as np
import pandas as pd
import xarray as xr

import src.Teleconnection.tools as tools
import src.Teleconnection.vertical_eof as vertical_eof


###### high level APIs #################
def season_eof(
    xarr: xr.DataArray,  # the geopotential data to be decomposed
    nmode: int = 2,  # the number of modes to generate
    window: int = 10,  # the rolling window if fixed_pattern = 'False'.
    fixed_pattern: str = "all",  # the fixed_pattern 'all','first','last','False'.
    independent: bool = True,  # the vertical eof strategy.
    standard: bool = True,  # standard pc or not.
):
    """high_level API for seasonal data eof analysis.

    **Arguments:**

        *xarr*: the DataArray to be decomposed. [time,ens,lat,lon,plev]
        *standard*: do standardization before the eof decompose or not.
        *independent*: all layers decompose independently or not.
        *rolling_eof*: whether to use rolling_eof or not.

    **Return:**

        EOF, PC and FRA
    """

    # passing parameters
    kwargs = {
        "nmode": nmode,  # for doeof
        "window": window,  # for rolling_eof
        "fixed_pattern": fixed_pattern,  # for rolling_eof
        "independent": independent,  # choose vetrical eof method.
        "standard": standard,
    }

    eof_result = vertical_eof.vertical_eof(xarr, **kwargs)

    return eof_result


def read_data(gph_dir, standardize=False):
    """
    read the gph data and do some pre-process.
    """
    # read MPI_onepct data
    try:
        zg_data = xr.open_dataset(gph_dir + "allens_zg.nc")
        zg_data = tools.split_ens(zg_data)
    except FileNotFoundError:
        zg_data = xr.open_mfdataset(
            gph_dir + "*.nc", combine="nested", concat_dim="ens", join="override"
        )
        zg_data = tools.split_ens(zg_data)
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
    if standardize:  # only standardize the data when decompose dependently.
        zg_trop = (zg_trop - zg_trop.mean(dim="time")) / zg_trop.std(dim="time")
    return zg_trop


def main():
    """
    for debug
    """
    # ex = xr.open_dataset("/work/mh0033/m300883/3rdPanel/data/sample.nc")
    # ex = ex.var156
    allens = xr.open_dataset(
        "/work/mh0033/m300883/transition/gr19/gphSeason/allens_season_time.nc"
    )
    # split ens
    splitens = tools.split_ens(allens)
    # demean ens-mean
    demean = splitens - splitens.mean(dim="ens")
    # select traposphere
    trop = demean.sel(plev=slice(85000, 100000)).isel(time=slice(0, 40))

    #     eof_sar,pc_sar,fra_sar = season_eof(ex,nmode=2,method ="rolling_eof",
    # window=10,fixed_pattern='all',return_full_eof= False,independent = True,standard=True)

    eof = season_eof(
        trop.var156,
        nmode=2,
        window=10,
        fixed_pattern="all",
        independent=True,
        standard=True,
    )


if __name__ == "__main__":
    main()
