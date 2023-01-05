import xarray as xr
import numpy as np
import src.composite.composite as composite


def Tel_field_composite(
    index: xr.DataArray,
    data: xr.DataArray,
    reduction: str = "mean",
    threshold: int = 2,
):
    """
    composite mean maps or counts of field in terms of teleconnection mode extremes.
    given the index and the data, find the coordinates (time and ens)
    of the extremes from the index, select the data with
    those coordinates, average/count the selected fields.
    **Arguments**
        *index* the index of NAO and EA
        *data* the original geopotential data.
        *reduction* mean or count
        *period* 'first10','last10','all'
    **Return**
        *compostie* the composite mean of the extreme cases.
    """

    # Select the same time period
    data = data.sel(time=index.time, method="nearest")

    # combine time and ens into one dim
    index = index.stack(com=("time", "ens"))
    data = data.stack(com=("time", "ens"))

    # since there is 'mode' dim in index, here groupby.
    tel_composite = index.groupby("mode").apply(
        composite.extreme_composite,
        data=data,
        reduction=reduction,
        dim="com",
        threshold=threshold,
    )

    return tel_composite


def field_composite(var, independent="dep", hlayer=100000):
    """
    function to get the first and last composite mean field
    """
    dataname = (
        "/work/mh0033/m300883/3rdPanel/data/influence/"
        + var
        + "/"
        + "onepct_1850-1999_ens_1-100."
        + var
        + ".nc"
    )

    indexname = (
        "/work/mh0033/m300883/3rdPanel/data/allPattern/"
        + independent
        + "_index_nonstd.nc"
    )

    changing_name = (
        "/work/mh0033/m300883/3rdPanel/data/changingPattern/"
        + independent
        + "_index_nonstd.nc"
    )

    # Data
    file = xr.open_dataset(dataname)
    fdata = file[var]
    # demean (ens-mean)
    demean = fdata - fdata.mean(dim="ens")

    # index
    all_all_dep = xr.open_dataset(indexname).pc
    changing_dep = xr.open_dataset(changing_name).pc
    all_all_dep = all_all_dep.transpose("time", "ens", "mode", "hlayers")

    mean_dep = all_all_dep.mean(dim="time")
    std_dep = all_all_dep.std(dim="time")
    dep_std = (changing_dep - mean_dep) / std_dep
    index = dep_std.sel(hlayers=hlayer)

    # change time
    index["time"] = index.time.dt.year
    demean["time"] = demean.time.dt.year

    ComFirst = composite(
        index=index, data=demean, dim="mode", reduction="mean", period="first10"
    )
    ComLast = composite(
        index=index, data=demean, dim="mode", reduction="mean", period="last10"
    )
    return ComFirst, ComLast, ComLast - ComFirst
