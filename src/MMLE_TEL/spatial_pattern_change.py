#%%
import xarray as xr
import numpy as np
import pandas as pd

#%%
import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools
import proplot as pplt
import cartopy.crs as ccrs


def read_gph_data(dir):
    data = xr.open_mfdataset(
        dir + "*.nc",
        combine="nested",
        concat_dim="ens",
        join="override",
    )
    data = data.rename({"plev": "plev"})
    try:
        data = data.zg
    except AttributeError:
        data = data.var156

    data = data.sel(plev=slice(100000, 20000))

    return data


# %%
def spatial_pattern_change(data, periods, names):
    """
    get the spatial pattern of the data in periods
    """
    if data.plev.size > 1:
        data = tools.standardize(data, dim="time")

    EOFs = []
    FRAs = []
    period_index = xr.IndexVariable(dims="period", data=periods)
    for period in periods:
        data_p = data.sel(time=period)
        EOF, FRA = vertical_spatial_pattern(data_p)
        EOFs.append(EOF)
        FRAs.append(FRA)
    EOFs = xr.concat(EOFs, dim=period_index)
    EOFs["period"] = ["0K", "2K", "4K"]

    FRAs = xr.concat(FRAs, dim=period_index)
    FRAs["period"] = names
    return EOFs, FRAs


def vertical_spatial_pattern(data, dim="plev"):
    """
    get the spatial pattern for all levels.
    """
    data = data.stack(com=("time", "ens"))
    eofs = data.groupby(dim, squeeze=True).apply(spatial_pattern, output="eof")
    fras = data.groupby(dim, squeeze=True).apply(spatial_pattern, output="fra")

    return eofs, fras


def spatial_pattern(data, output="eof"):
    """
    get the spatial pattern of single data
    """
    eof, _, fra = ssp.doeof(data, nmode=2, dim="com", standard=True)
    if output == "eof":
        out = eof
    elif output == "fra":
        out = fra
    return out


# %%


def spatial_stat(eof, mode, dim="lon"):
    """
    zonally mean ('lat' as the final dim) or meridional mean ('lon' as the final dim)
    **Arguments**
        *eof* the eof
        *mode* the mode
        *dim* the dim that left
    """
    if dim == "lon":
        bins = np.arange(-90, 41, 5)  # lon bins
        labels = np.arange(-90, 36, 5)  # lon label
        average_dim = "lat"
    elif dim == "lat":
        bins = np.arange(20, 81, 4)
        labels = np.arange(22, 81, 4)
        average_dim = "lon"

    data = eof.sel(mode=mode)
    data_height = data.groupby_bins(dim, bins=bins, labels=labels).mean(
        dim=average_dim, skipna=True
    )
    return data_height


# %%


# PLOT maps
def spatial_pattern_maps(eofs, fras, levels=np.arange(-1.0, 1.1, 0.2)):
    """
    rows as modes
    cols in different periods
    """
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
        leftlabels=("NAO", "EA"),
        suptitle=f"spatial patterns on 500hPa",
    )

    for r, mode in enumerate(eofs.mode):
        for c, period in enumerate(eofs.period):

            eof = eofs.sel(mode=mode, period=period)
            fra = fras.sel(mode=mode, period=period)
            map = axes[r, c].contourf(
                eof,
                x="lon",
                y="lat",
                levels=levels,
                extend="both",
                transform=ccrs.PlateCarree(),
                cmap="RdBu_r",
            )
            axes[r, c].set_title(str(period.values) + f"({fra.values:.0%})")

    fig.colorbar(map, loc="r", pad=3, title="gph/std")


# %%


def spatial_pattern_profile(eofs, levels=np.arange(-2.0, 2.1, 0.4)):
    """
    rows mode
    cols periods
    """
    eofs["plev"] = eofs["plev"] / 100

    lon_NAO = spatial_stat(eofs, mode="NAO", dim="lon")

    lat_EA = spatial_stat(eofs, mode="EA", dim="lat")

    fig, axes = pplt.subplots(
        space=0,
        refwidth="25em",
        wspace=3,
        hspace=3,
        nrows=2,
        ncols=3,
        share=False,
    )
    axes.format(
        latlines=20,
        lonlines=30,
        toplabels=("0K", "2K", "4K"),
        leftlabels=("NAO", "EA"),
        suptitle=f"spatial change profile",
        ylim=(1000, 200),
        ylabel="gph/hpa",
    )

    xs = ["lon", "lat"]

    for r, profile in enumerate([lon_NAO, lat_EA]):
        for c, period in enumerate(eofs.period):
            prof = profile.sel(period=period)
            vertmap = axes[r, c].contourf(
                prof, x=xs[r], y="plev", levels=levels, extend="both", cmap="RdBu_r"
            )
    fig.colorbar(vertmap, loc="r", pad=3, title="gph/std")


# %%
# given a time slice called period, for example ['1851-12-31','1860-12-31'],find the middle time
def get_middle_time(period):
    """
    get the middle time of a period
    """
    start = np.array(period.start).astype("datetime64[Y]")
    end = np.array(period.stop).astype("datetime64[Y]")
    middle = start + (end - start) / 2
    return middle


# plot the spatial pattern changes under 0K,2K,4K warming
def spatial_pattern_change_decade(
    periods, patterns, fras, plev=50000, levels=np.arange(-1.0, 1.1, 0.2)
):
    # data preapration
    patterns = patterns.sel(plev=plev)
    fras = fras.sel(plev=plev)
    middle_years = [get_middle_time(period) for period in periods]
    patterns_warming = patterns.sel(decade=middle_years, method="nearest")

    # plot
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
        leftlabels=("NAO", "EA"),
        suptitle=f"spatial patterns on {plev/100:.0%}hPa",
    )

    period_names = ["0K", "2K", "4K"]
    for r, mode in enumerate(patterns_warming.mode):
        for c, period in enumerate(patterns_warming.decade):

            eof = patterns_warming.sel(mode=mode, decade=period)
            fra = fras.sel(mode=mode, decade=period)
            map = axes[r, c].contourf(
                eof,
                x="lon",
                y="lat",
                levels=levels,
                extend="both",
                transform=ccrs.PlateCarree(),
                cmap="RdBu_r",
            )
            axes[r, c].set_title(
                period_names[c]
                + "-"
                + str(period.dt.year.values)
                + f"({fra.values:.0%})"
            )

    fig.colorbar(map, loc="r", pad=3, title="gph/std")
    return fig
