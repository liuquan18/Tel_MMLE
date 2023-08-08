#%%
# import xarray, cartopy, matplotlib
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import proplot as pplt
import seaborn as sns
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker

import src.warming_stage.warming_stage as warming_stage

# no warnings
import warnings

warnings.filterwarnings("ignore")


# Fig 1 spatial patterns and statistics of the pcs
# rows for 'NAO' and 'EA'
# colums for 'spatial map','pc hist','violion vertical profile'.


#%% for one model spatial pattern and index distritbuion together 
def stat_overview(
    first_pc,
    last_pc,
    first_fra=None,
    last_fra=None,
    first_eof=None,
    last_eof=None,
    levels=np.arange(-2, 2.1, 0.4),
):

    # plot
    fig = pplt.figure(space=0, refwidth="25em", wspace=3, hspace=3)
    fig.format(
        abc=True,
        abcloc="ul",
        abcstyle="a",
        title="spatial patterns and statistics of the pcs",
        leftlabels=("NAO", "EA"),
    )

    gs = pplt.GridSpec(
        ncols=2,
        nrows=2,
        wspace=2,
        wratios=(1, 1),
    )
    modes = ["NAO", "EA"]

    for i, mode in enumerate(modes):
        # data preparation

        ## eof as xr.DataArray
        first_eof_mode = first_eof.sel(mode=mode).squeeze()
        first_fra_mode = first_fra.sel(mode=mode).squeeze()

        try:
            last_eof_mode = last_eof.sel(mode=mode).squeeze()
            last_fra_mode = last_fra.sel(mode=mode).squeeze()
        except AttributeError:
            pass

        spatial_ax = fig.add_subplot(gs[i, 0],proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})

        # plot spatial map at 500hPa
        spatial_ax, fmap, lmap = spatial_pattern_plot(
            spatial_ax,
            first_eof_mode,
            first_fra_mode,
            last_eof_mode,
            last_fra_mode,
            levels,
        )

        # plot pc hist
        # split the pc into two parts
        first_pc_mode = first_pc.sel(mode=mode)
        last_pc_mode = last_pc.sel(mode=mode)
        hist_ax = fig.add_subplot(gs[i, 1])
        hist_ax = index_distribution_plot(hist_ax, first_pc_mode, last_pc_mode)

        # add legend
        f_patch = mpatches.Patch(color="#1f77b4", label="first10")
        l_patch = mpatches.Patch(color="#ff7f0e", label="last10")

        if i == 1:
            spatial_ax.colorbar(fmap, loc="b", title="std", ticks=0.3, pad=2)
            hist_ax.legend(
                handles=[f_patch, l_patch], loc="b", title="periods", frameon=False
            )
    return fig



#%% for all MMLEAS, plot the spatial pattern and index distribution
def spatial_index_MMLEA(eof_firsts, eof_lasts):
    models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    models_legend = [
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    ]

    fig = pplt.figure(space=0, refwidth="50em", wspace=3, hspace=3)

    gs = pplt.GridSpec(
    ncols=5,
    nrows=2,
    hratios=(
        1.4,
        1,
    ),
)

# rows for spatial pattern and index, columns for different models
    for c, model in enumerate(models):
    # plot spatial pattern
        spatial_ax = fig.add_subplot(gs[0, c],proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
        spatial_ax, fmap, lmap = spatial_pattern_plot(
        spatial_ax,
        eof_firsts[c].eof,
        eof_firsts[c].fra,
        eof_lasts[c].eof,
        eof_lasts[c].fra,
    )

    # plot index distribution
        hist_ax = fig.add_subplot(gs[1, c])
        hist_ax,hist = index_distribution_plot(
        hist_ax,
        eof_firsts[c].pc,
        eof_lasts[c].pc,
    )

    # add legend
        f_patch = mpatches.Patch(color="#1f77b4", label="first10")
        l_patch = mpatches.Patch(color="#ff7f0e", label="last10")

    spatial_ax.colorbar(fmap, loc="r", title="std", pad=2)
    hist_ax.legend(
    handles=[f_patch, l_patch], loc="r", title="periods", frameon=False,
    ncols = 1,
)
    fig.format(
    suptitle="spatial pattern and distribution of summer NAO at 500hpa",
    sharex=False,
    sharey=False,
    bottomlabels = models_legend,
)
    return fig

# %% tool function for plotting the spatial pattern 
def spatial_pattern_plot(
    spatial_ax,
    first_eof,
    first_fra,
    last_eof=None,
    last_fra=None,
    levels=np.arange(-2, 2.1, 0.4),
    title="",
):

    fmap = spatial_ax.contourf(
        first_eof.lon,
        first_eof.lat,
        first_eof,
        levels=levels,
        extend="both",
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
    )

    lmap = None
    if last_eof is not None:
        lmap = spatial_ax.contour(
            last_eof.lon,
            last_eof.lat,
            last_eof,
            colors="gray8",
            nozero=True,
            labels=True,
            levels=np.delete(levels, int((len(levels) - 1) / 2)),
            labels_kw={"weight": "bold"},
            add_colorbar=False,
            lw = 1.0,
        )

        spatial_ax.format(
            lonlines=20,
            latlines=30,
            coast=True,
            coastlinewidth=0.5,
            coastcolor="charcoal",
            title=f"({first_fra:.0%}" + f"->{last_fra:.0%})",
        )
    else:
        spatial_ax.format(
            lonlines=20,
            latlines=30,
            coast=True,
            coastlinewidth=0.5,
            coastcolor="charcoal",
            title=f"{title} ({first_fra:.0%})",
        )
    return spatial_ax, fmap, lmap


#%% tool function to plot the index distribution
def index_distribution_plot(hist_ax, first_pc, last_pc):
    df = index_to_df(first_pc, last_pc)
    hist = sns.histplot(
        data=df,
        x="pc",
        hue="periods",
        hue_order=["first", "last"],
        palette=["#1f77b4", "#ff7f0e"],
        multiple="dodge",
        shrink=0.6,
        bins=np.arange(-4, 4.1, 0.5),
        legend=False,
        ax=hist_ax,
        stat= "density",
    )

    hist_ax.format(grid=False, 
                   yminorticks="null", 
                   xminorticks="null",
                   xticks = [-3,-1.5,0,1.5,3])
    hist_ax.spines["right"].set_visible(False)
    hist_ax.spines["top"].set_visible(False)
    return hist_ax,hist


#%% tool function to plot the lines of obs and mmlea (envelop or not)
def envelop_obs_mmlea(ax, obs, mmlea):

    max_mmlea = mmlea.max(dim="ens")
    min_mmlea = mmlea.min(dim="ens")
    # shading between the max and min
    ax.fill_between(
        np.arange(max_mmlea.time.size),
        max_mmlea.values,
        min_mmlea.values,
        color="grey7",
        alpha=0.3,
        linewidth=0,
    )

    ax.plot(obs.values, label="obs", color="orange", linewidth=2)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(30))
    # ticklabels = [str(year) for year in np.insert(obs.time[0::30].dt.year.values, 0, 0)]
    # ax.xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_year_summer))
    ax.format(
        ylim=(-3.5, 3.5),
        yminorticks="null",
        xminorticks="null",
        grid=False,
    )
    return ax

def format_year_summer(value, tick_number):
    value_year = int(value/3) + 1940
    return f"JJA \n {str(int(value_year))}"

def index_to_df(first_pc, last_pc):
    periods = xr.IndexVariable("periods", ["first", "last"])
    both = xr.concat([first_pc, last_pc], dim=periods)
    both = both.to_dataframe().reset_index()
    both = both[["periods", "pc"]]
    return both

#%% tool function to plot the box plot of obs and mmlea
def obs_mmlea_box_plot( box_ax,EOFs):
    models_legend = [
        "ERA5 (obs)",
        "MPI-GE (100)",
        "CanESM2 (50)",
        "CESM1 (40)",
        "MK3.6 (30)",
        "GFDL-CM3 (20)",
    ]
    models = ["ERA5", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    keys = EOFs.keys()
    colors = ['w', "C1", "tab:purple", "tab:blue", "tab:green", "C4"]
    fills = [False, True, True, True, True, True]
    bps = []
    for i, key in enumerate(models):
        bps.append(box_ax.boxplot(
        i + 1,
        EOFs[key].pc.values.reshape(-1),
        flierprops={"markerfacecolor": "grey", "marker": "o", "markersize": 0.8},
        widths=0.5,
        fc = colors[i],
        fill = fills[i],
        label = models_legend[i],
        alpha = 0.8,
    ))
        box_ax.set_ylim(-3.5, 3.5)

    box_ax.format(
    xminorticks="null",
    yminorticks="null",
    grid=False,
)

    box_ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    box_ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    return box_ax

def format_func(value, tick_number):
    if value == 1:
        return "ERA5 \n (obs)"
    elif value == 2:
        return "MPI-GE \n (100)"
    elif value == 3:
        return "CanESM2 \n (50)"
    elif value == 4:
        return "CESM1-CAM5 \n (40)"
    elif value == 5:
        return "MK3.6 \n (30)"
    elif value == 6:
        return "GFDL-CM3 \n (20)"
    else:
        return ""
    
# %%
