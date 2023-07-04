import xarray as xr
import numpy as np
import src.composite.composite as composite_analysis
import cartopy.crs as ccrs
import proplot as pplt
import src.plots.utils as utils

def Tel_field_composite(
    index: xr.DataArray,
    data: xr.DataArray,
    threshold: float = 2,
    reduction = "mean",
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
    index_c = index.copy() # make a copy of the original data
    data_c = data.copy()

    index_c['time'] = index_c.time.dt.year
    data_c['time'] = data_c.time.dt.year
    data_c = data_c.sel(time=index_c.time, method="nearest")

    # combine time and ens into one dim
    index_c = index_c.stack(com=("time", "ens"))
    data_c = data_c.stack(com=("time", "ens"))

    # since there is 'mode' dim in index, here groupby.
    tel_composite = index_c.groupby("mode").apply(
        composite_analysis.extreme_composite,
        data=data_c,
        dim="com",
        threshold=threshold,
        reduction=reduction,
    )

    return tel_composite


def composite_plot( first, last, mode,level_bound = 2):
    if mode == "NAO":
        bound_l = -1*level_bound -1
        bound_u = level_bound +1 + 0.1
    elif mode == "EA":
        bound_l = -1*level_bound
        bound_u = level_bound + 0.1
    levels = np.arange(bound_l,bound_u,0.5)

    first = utils.erase_white_line(first)
    last = utils.erase_white_line(last)

    data_all = [
        first.sel(mode=mode),
        last.sel(mode=mode),
        last.sel(mode=mode) - first.sel(mode=mode),
    ]
    extr_type = ["pos", "neg"]

    fig, axes = pplt.subplots(
        space=0,
        refwidth="25em",
        wspace=3,
        hspace=3,
        proj="ortho",
        proj_kw=({"lon_0": -20, "lat_0": 60}),
        nrows=2,
        ncols=3,
    )
    axes.format(
        latlines=20,
        lonlines=30,
        coast=True,
        coastlinewidth=0.5,
        coastcolor="gray7",
        toplabels=["first10", "last10", "last10 - first10"],
        leftlabels=("pos", "neg"),
        suptitle=f"Change in influence of extreme {mode} on surface temperature",
    )

    extr_types = ["pos", "neg"]
    for i, extr_type in enumerate(extr_types):
        for j, data in enumerate(data_all):  # one row

            first_m = axes[i, j].contourf(
                data.sel(extr_type=extr_type),
                x="lon",
                y="lat",
                levels=levels,
                extend="both",
                transform=ccrs.PlateCarree(),
                cmap="RdBu_r",
            )

    fig.colorbar(first_m, loc="r", pad=3, title=f"tsurf/K")



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
    all_all_dep = all_all_dep.transpose("time", "ens", "mode", "plev")

    mean_dep = all_all_dep.mean(dim="time")
    std_dep = all_all_dep.std(dim="time")
    dep_std = (changing_dep - mean_dep) / std_dep
    index = dep_std.sel(plev=hlayer)

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
