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
    both["plev"] = (both["plev"] / 100).astype(int)
    return both


def stat_overview(eof_result, plev=50000):

    eof = eof_result.eof
    pc = eof_result.pc
    fra = eof_result.fra

    # split eof_result into first 10 and last 10 years
    first_eof, last_eof, first_pc, last_pc, first_fra, last_fra = split_period(
        eof, pc, fra
    )

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
    levels = np.arange(-0.8, 0.9, 0.2)

    for i, mode in enumerate(modes):
        # data preparation
        ## eof as xr.DataArray
        first_eof_500 = first_eof.sel(mode=mode, plev=50000).squeeze()
        last_eof_500 = last_eof.sel(mode=mode, plev=50000).squeeze()

        ## fras selecte the plev
        try:
            first_fra_500 = first_fra.sel(mode=mode, plev=50000).squeeze()
            last_fra_500 = last_fra.sel(mode=mode, plev=50000).squeeze()
        except KeyError:
            first_fra_500 = first_fra.sel(mode=mode).squeeze()
            last_fra_500 = last_fra.sel(mode=mode).squeeze()

        ## pc to dataframe
        df = to_dataframe(first_pc, last_pc, mode)
        df_500 = df[df["plev"] == 500]

        # plot spatial map at 500hPa
        spatial_ax = fig.add_subplot(
            gs[i, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60}
        )

        fmap = first_eof_500.plot.contourf(
            ax=spatial_ax, levels=levels, extend="both", add_colorbar=False
        )

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
            title=mode
            + f"({first_fra_500:.0%}"
            + f"->{last_fra_500:.0%})",
        )

        # plot pc hist
        hist_ax = fig.add_subplot(gs[i, 1])

        hist = sns.histplot(
            data=df_500,
            x=mode,
            hue="compare",
            hue_order=["first10", "last10"],
            palette=['#1f77b4', '#ff7f0e'],
            multiple="dodge",
            shrink=0.6,
            bins=np.arange(-4, 4.1, 0.5),
            legend=False,
            ax=hist_ax,
        )

        hist_ax.format(grid=False, yminorticks="null", xminorticks="null", title=mode)
        hist_ax.spines["right"].set_visible(False)
        hist_ax.spines["top"].set_visible(False)

        # plot violin
        violin_ax = fig.add_subplot(gs[i, 2])
        g = sns.violinplot(
            data=df,
            y="plev",
            x=modes[i],
            hue="compare",
            kind="violin",
            palette=['#1f77b4', '#ff7f0e'],
            orient="h",
            ax=violin_ax,
            split=False,
            dodge=True,
            linewidth=1,
            alpha=0.3,
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

        # add legend
        f_patch = mpatches.Patch(color="#1f77b4", label="first10")
        l_patch = mpatches.Patch(color="#ff7f0e", label="last10")

        if i == 1:
            spatial_ax.colorbar(fmap, loc="b", title="std", ticks=0.2, pad=2)
            hist_ax.legend(
                handles=[f_patch, l_patch], loc="b", title="periods", frameon=False
            )
            violin_ax.legend(loc="b", ncols=2, title="periods", frameon=False)

    return fig


# %%
