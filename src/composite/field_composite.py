#%%
import xarray as xr
import numpy as np
import src.composite.composite as composite_analysis
import cartopy.crs as ccrs
import proplot as pplt
import src.plots.utils as utils

def Tel_field_composite(
    index: xr.DataArray,
    data: xr.DataArray,
    threshold: float = 1.5,
    reduction = "mean",
    **kwargs
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

    # select the same time period
    data_c = data_c.sortby("time")
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
        **kwargs
    )

    return tel_composite


def composite_plot( first, last, mode,level_bound = None,levels = None):
    if levels is None:
        if mode == "NAO":
            bound_l = -1*level_bound -1
            bound_u = level_bound +1 + 0.1
            interval = level_bound/4
        elif mode == "EA":
            bound_l = -1*level_bound
            bound_u = level_bound + 0.1
            interval = level_bound/5
        levels = np.arange(bound_l,bound_u,interval)
    else:
        levels = levels

    first = utils.erase_white_line(first)
    last = utils.erase_white_line(last)

    data_all = [
        first.sel(mode=mode),
        last.sel(mode=mode),
        last.sel(mode=mode) - first.sel(mode=mode),
    ]
    extr_type = ["pos", "neg"]

    params = {
        "ytick.color": "w",
        "xtick.color": "w",
        "axes.labelcolor": "w",
        "axes.edgecolor": "w",
        "tick.labelcolor": "w",
        "text.color": "w",
        "fontsize": 25,
        "figure.facecolor": "black",
        "axes.facecolor": "black",
        "grid.linewidth": 1.5,
        "grid.color": "black",
    }

    fig, axes = pplt.subplots(
        space=0,
        refwidth="35em",
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
        color = 'grey7',
        coast=True,
        coastlinewidth=1,
        coastcolor="grey9",
        toplabels=["first10", "last10", "last10 - first10"],
        toplabelcolor="w",
        toplabels_kw = {"fontsize": 50},
        leftlabels=("pos", "neg"),
        leftlabelcolor="w",
        leftlabels_kw = {"fontsize": 50},
        suptitle=f"Change in influence of extreme {mode} on surface temperature",
        # set the fontsize of labels to 25        
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
            axes[i,j].grid(color = 'grey7',linewidth = 1.5)
    fig.colorbar(first_m, loc="r", pad=3, title=f"tsurf/K")


# %%
def composite_plot_MMLEA(firsts,lasts,mode='NAO',level_bound = None,levels = np.arange(-1,1.1,0.2),extr_type = 'pos'):
    models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    models_legend = [
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    ]

    fig, axes = pplt.subplots(
        space=0,
        refwidth="50em",
        wspace=3,
        hspace=3,
        proj="ortho",
        proj_kw=({"lon_0": -20, "lat_0": 60}),
        nrows=3,
        ncols=5,
    )
    axes.format(
        latlines=20,
        lonlines=30,
        color = 'grey7',
        coast=True,
        coastlinewidth=1,
        coastcolor="grey9",
        leftlabels=["first10", "last10", "last10 - first10"],
        toplabels=models_legend,
        toplabels_kw = {"fontsize": 50},
        leftlabels_kw = {"fontsize": 50},
        suptitle=f"Change in influence of extreme {mode} on surface temperature",
        # set the fontsize of labels to 25        
    )


    for i, model in enumerate(models): # cols for different models
        first = firsts[model].sel(mode=mode,extr_type = extr_type)
        last = lasts[model].sel(mode=mode,extr_type = extr_type)
        first = utils.erase_white_line(first)
        last = utils.erase_white_line(last)

        data_all = [
            first,
            last,
            last - first
        ]

        for j, data in enumerate(data_all):  # row for different data
            first_m = axes[j, i].contourf(
                data,
                x="lon",
                y="lat",
                levels=levels,
                extend="both",
                transform=ccrs.PlateCarree(),
                cmap="RdBu_r",
            )
            axes[j,i].grid(color = 'grey7',linewidth = 1.5)
    fig.colorbar(first_m, loc="r", pad=3, title=f"tsurf/K")



def plot_composite_single_ext(COMPOSITEs, models, axes,extr_type = 'pos'):
    for i, model in enumerate(models): # cols for different models
        first = COMPOSITEs[model][0].sel(mode='NAO',extr_type = extr_type)
        last = COMPOSITEs[model][-1].sel(mode='NAO',extr_type = extr_type)
        first = utils.erase_white_line(first)
        last = utils.erase_white_line(last)

        data_all = [
        first,
        last,
        last - first
    ]

        maps = []
        for j, data in enumerate(data_all):  # row for different data
            map = axes[j, i].contourf(
            data,
            x="lon",
            y="lat",
            levels = np.arange(-1.5,1.6,0.3),
            extend="both",
            transform=ccrs.PlateCarree(),
            cmap="RdBu_r",
        )
            maps.append(map)
            axes[j,i].grid(color = 'grey7',linewidth = 0.5)
    return axes, maps
