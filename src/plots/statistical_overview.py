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


# Fig 1 spatial patterns and statistics of the pcs
# rows for 'NAO' and 'EA'
# colums for 'spatial map','pc hist','violion vertical profile'.


def split_period(eof, pc, fra):
    # split eof_result into first 10 and last 10 years
    periods_pc, periods = warming_stage.split_period(pc, compare="CO2")
    first_pc, last_pc = periods_pc[0], periods_pc[1]
    first_eof = eof.sel(decade=periods[0])
    last_eof = eof.sel(decade=periods[1])
    first_fra = fra.sel(decade=periods[0])
    last_fra = fra.sel(decade=periods[1])
    return first_eof, last_eof, first_pc, last_pc, first_fra, last_fra


def to_dataframe(first_pc, last_pc, mode):
    first = first_pc.sel(mode=mode).to_dataframe(mode)
    last = last_pc.sel(mode=mode).to_dataframe(mode)
    both = pd.concat([first, last], axis=0)
    both = both.reset_index()
    both = both[["plev", "compare", mode]]
    both['plev'] = (both['plev']/100).astype(int)
    return both


def stat_overview(eof_result, plev=50000):

    # select plev
    eof = eof_result.eof
    pc = eof_result.pc
    fra = eof_result.fra

    # split eof_result into first 10 and last 10 years
    first_eof, last_eof, first_pc, last_pc, first_fra, last_fra = split_period(
        eof, pc, fra
    )

    # pc to dataframe
    coords = xr.IndexVariable(dims="periods", data=["first10", "last10"])
    index_all_periods = xr.concat([first_pc, last_pc], dim=coords)
    index_all_periods = index_all_periods.to_dataframe().reset_index()

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
        ncols=3,
        nrows=2,
        wspace=2,
        wratios=(1, 1, 1),
    )
    modes = ["NAO", "EA"]
    levels = np.arange(-1, 1.1, 0.2)

    for i, mode in enumerate(modes):
        # plot spatial map at 500hPa
        spatial_ax = fig.add_subplot(
            gs[i, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60}
        )

        first_data = first_eof.sel(mode=mode, plev=50000).squeeze()
        last_data = last_eof.sel(mode=mode, plev=50000).squeeze()

        fmap = first_data.plot.contourf(
            ax=spatial_ax, levels=levels, extend="both", add_colorbar=False
        )

        lmap = last_data.plot.contour(
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
            title=mode
            + f"({first_fra.sel(mode = mode,plev = 50000).squeeze():.0%}"
            + f"->{last_fra.sel(mode = mode,plev = 50000).squeeze():.0%})",
        )

        # plot pc hist
        hist_ax = fig.add_subplot(gs[i, 1])

        # if i == 1, add legend
        if i == 1:
            legend = True
        else:
            legend = False

        hist = sns.histplot(
            data=index_all_periods[index_all_periods["mode"] == mode],
            x="pc",
            hue="periods",
            hue_order=["first10", "last10"],
            multiple="dodge",
            shrink=1,
            bins=np.arange(-4, 4.1, 0.5),
            legend = legend,
        )

        # legend
        fpatch = mpatches.Patch(color="grey", label="first10 years")
        lpatch = mpatches.Patch(color="black", label="first10 years")
        handles = [fpatch, lpatch]

        hist_ax.format(grid=False, yminorticks="null", xminorticks="null", title=mode)
        hist_ax.spines["right"].set_visible(False)
        hist_ax.spines["top"].set_visible(False)

        # plot violin
        violin_ax = fig.add_subplot(gs[i, 2])
        df = to_dataframe(first_pc, last_pc, mode)
        g = sns.violinplot(
            data=df,
            y="plev",
            x=modes[i],
            hue="compare",
            kind="violin",
            palette="pastel",
            orient="h",
            ax=violin_ax,
            split=False,
            dodge=True,
            linewidth=1,
        )
        g.axes.legend().remove()
        violin_ax.format(
            grid=False,
            yminorticks="null",
            xminorticks="null",
            title=mode,
            xlabel="std",
            ylabel="gph/hPa",
            xlim=(-5, 5),
        )
        violin_ax.spines["left"].set_visible(False)
        violin_ax.spines["right"].set_visible(False)
        violin_ax.spines["top"].set_visible(False)

        if i == 1:
            spatial_ax.colorbar(fmap, loc="b", title="std", ticks=0.2, pad=2)
            sns.move_legend(hist_ax, "b", title="periods", ncol=2)
            violin_ax.legend(loc="b", ncols=2, title="periods")


    return fig


# %%
