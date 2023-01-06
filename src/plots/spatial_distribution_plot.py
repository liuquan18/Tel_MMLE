import proplot as pplt
import seaborn as sns
import numpy as np
import matplotlib as plt

# plot
def spatialMap_violin(eof_500hpa, index_500hpa, fra_500hpa):
    """
    sptail maps and violin plots of index
    """
    fig = pplt.figure(space=0, refwidth="25em", wspace=3, hspace=3)

    fig.format(
        suptitle="spatial pattern and distribution of NAO and EA at 500hpa",
    )

    gs = pplt.GridSpec(
        ncols=2,
        nrows=2,
        hratios=(
            1.4,
            1,
        ),
    )
    modes = ["NAO", "EA"]

    for i, mode in enumerate(modes):
        ax1 = fig.subplot(gs[0, i], proj="ortho", proj_kw=({"lon_0": -20, "lat_0": 60}))

        ax1.format(
            latlines=20,
            lonlines=30,
            coast=True,
            coastlinewidth=0.5,
            coastcolor="charcoal",
            title=modes[i] + f" ({fra_500hpa.sel(mode = mode):.0%})",
        )

        map = ax1.contourf(
            eof_500hpa.sel(mode=mode), levels=np.arange(-1, 1.1, 0.1), extend="both"
        )
        if i == 1:
            ax1.colorbar(map, loc="r", title="std", ticks=0.2, pad=2)

        ax2 = fig.subplot(gs[1, i])
        violin = sns.violinplot(
            data=index_500hpa[index_500hpa["mode"] == mode],
            y="pc",
            x="periods",
            ax=ax2,
            orient="v",
        )

        ax2.format(
            ylim=(-4.5, 4.5),
            ytickminor=False,
            xlocator=(0, 1),
        )

        ax2.spines.right.set_visible(False)
        ax2.spines.top.set_visible(False)
    return fig


# %%
# plot
def spatialMap_hist(eof_500hpa, index_500hpa, fra_500hpa):
    """
    spatial maps and hist of index
    """
    fig = pplt.figure(space=0, refwidth="25em", wspace=3, hspace=3)

    fig.format(
        suptitle="spatial pattern and distribution of NAO and EA at 500hpa",
    )

    gs = pplt.GridSpec(
        ncols=2,
        nrows=2,
        hratios=(
            1.4,
            1,
        ),
    )
    modes = ["NAO", "EA"]

    for i, mode in enumerate(modes):
        ax1 = fig.subplot(gs[0, i], proj="ortho", proj_kw=({"lon_0": -20, "lat_0": 60}))

        ax1.format(
            latlines=20,
            lonlines=30,
            coast=True,
            coastlinewidth=0.5,
            coastcolor="charcoal",
            title=modes[i] + f" ({fra_500hpa.sel(mode = mode):.0%})",
        )

        map = ax1.contourf(
            eof_500hpa.sel(mode=mode), levels=np.arange(-1, 1.1, 0.1), extend="both"
        )
        if i == 1:
            ax1.colorbar(map, loc="r", title="std", ticks=0.2, pad=2)

        ax2 = fig.subplot(gs[1, i])
        violin = sns.histplot(
            data=index_500hpa[index_500hpa["mode"] == mode],
            x="pc",
            hue="periods",
            multiple="dodge",
            shrink=1,
            bins=np.arange(-4, 4.1, 0.5),
        )
        violin.legend(ncol=1, labels=["last10", "first10"], title="period")

        ax2.format(
            xtickminor=False,
            ytickminor=False,
            grid=False,
        )
        ax2.spines.right.set_visible(False)
        ax2.spines.top.set_visible(False)
    return fig


# %%
