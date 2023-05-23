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


def split_pc(first_pc, last_pc, mode):
    # split the pc into two parts
    periods = xr.IndexVariable("periods", ["first", "last"])

    first_pc = first_pc.sel(mode=mode)
    last_pc = last_pc.sel(mode=mode)

    both = xr.concat([first_pc, last_pc], dim=periods)
    both = both.to_dataframe().reset_index()
    both = both[["periods", "pc"]]
    return both


def stat_overview(
    first_eof = None,
    last_eof = None,
    first_pc,
    last_pc,
    first_fra,
    last_fra,
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
        first_eof_500 = first_eof.sel(mode=mode).squeeze()
        last_eof_500 = last_eof.sel(mode=mode).squeeze()

        first_fra_500 = first_fra.sel(mode=mode).squeeze()
        last_fra_500 = last_fra.sel(mode=mode).squeeze()

        ## pc to dataframe
        df_500 = split_pc(first_pc, last_pc, mode)

        # plot spatial map at 500hPa
        spatial_ax = fig.add_subplot(
            gs[i, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60}
        )

        fmap = first_eof_500.plot.contourf(
            ax=spatial_ax, levels=levels, extend="both", add_colorbar=False
        )
        if last_eof is not None:
            lmap = last_eof_500.plot.contour(
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
            title=mode + f"({first_fra_500:.0%}" + f"->{last_fra_500:.0%})",
        )

        # plot pc hist
        hist_ax = fig.add_subplot(gs[i, 1])

        hist = sns.histplot(
            data=df_500,
            x="pc",
            hue="periods",
            hue_order=["first", "last"],
            palette=["#1f77b4", "#ff7f0e"],
            multiple="dodge",
            shrink=0.6,
            bins=np.arange(-4, 4.1, 0.5),
            legend=False,
            ax=hist_ax,
        )

        hist_ax.format(grid=False, yminorticks="null", xminorticks="null", title=mode)
        hist_ax.spines["right"].set_visible(False)
        hist_ax.spines["top"].set_visible(False)

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
