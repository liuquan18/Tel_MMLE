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

import src.warming_stage.warming_stage as warming_stage

# no warnings
import warnings

warnings.filterwarnings("ignore")


# Fig 1 spatial patterns and statistics of the pcs
# rows for 'NAO' and 'EA'
# colums for 'spatial map','pc hist','violion vertical profile'.


def index_to_df(first_pc, last_pc):
    periods = xr.IndexVariable("periods", ["first", "last"])
    both = xr.concat([first_pc, last_pc], dim=periods)
    both = both.to_dataframe().reset_index()
    both = both[["periods", "pc"]]
    return both


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

        # plot spatial map at 500hPa
        spatial_ax, fmap, lmap = spatial_pattern_plot(
            fig,
            gs[i, 0],
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
        hist_ax = index_distribution_plot(fig, gs[i,1], first_pc_mode, last_pc_mode)

        # add legend
        f_patch = mpatches.Patch(color="#1f77b4", label="first10")
        l_patch = mpatches.Patch(color="#ff7f0e", label="last10")

        if i == 1:
            spatial_ax.colorbar(fmap, loc="b", title="std", ticks=0.3, pad=2)
            hist_ax.legend(
                handles=[f_patch, l_patch], loc="b", title="periods", frameon=False
            )
    return fig


# %%
def spatial_pattern_plot(
    fig,
    gs,
    first_eof,
    first_fra,
    last_eof=None,
    last_fra=None,
    levels=np.arange(-2, 2.1, 0.4),
):
    # plot spatial map at 500hPa
    spatial_ax = fig.add_subplot(gs, proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})

    fmap = first_eof.plot.contourf(
        ax=spatial_ax, levels=levels, extend="both", add_colorbar=False
    )

    lmap = None
    if last_eof is not None:
        lmap = last_eof.plot.contour(
            ax=spatial_ax,
            colors="gray8",
            nozero=True,
            labels=True,
            levels=np.delete(levels, int((len(levels) - 1) / 2)),
            labels_kw={"weight": "bold"},
            add_colorbar=False,
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
            title=f"({first_fra_500:.0%})",
        )
    return spatial_ax, fmap, lmap



def index_distribution_plot(fig, gs, first_pc, last_pc):
    hist_ax = fig.add_subplot(gs)
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
        spatial_ax, fmap, lmap = spatial_pattern_plot(
        fig,
        gs[0, c],
        eof_firsts[c].eof,
        eof_firsts[c].fra,
        eof_lasts[c].eof,
        eof_lasts[c].fra,
    )

    # plot index distribution
        hist_ax,hist = index_distribution_plot(
        fig,
        gs[1, c],
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