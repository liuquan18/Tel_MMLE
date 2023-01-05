import xarray as xr
import pandas as pd
import proplot as pplt
import numpy as np
import seaborn as sns


def to_dataframe(period, mode, name):
    """
    change single xarray as pd.DataFrame
    """
    period = period.sel(mode=mode).drop("mode")
    period = period.stack(com=("time", "ens"))
    period["com"] = np.arange(period.com.size)
    period = period.to_dataframe(name)
    return period


def xr2df(first10_all_period, last10_all_period, mode):
    """
    prepare dataframe for seaborn.
    """
    first = to_dataframe(first10_all_period, mode=mode, name=mode)
    last = to_dataframe(last10_all_period, mode=mode, name=mode)
    periods = pd.concat(
        [first, last], keys=["first10", "last10"], names=["period", "hlayers", "com"]
    )
    periods = periods.reset_index().set_index("com")
    periods["hlayers"] = periods["hlayers"] / 100
    return periods


def plot_vilion(first, last, std_type, split=False):
    fig = pplt.figure(space=0, refwidth="25em")
    axes = fig.subplots(nrows=1, ncols=2)
    axes.format(
        abc="a",
        abcloc="ul",
        xlocator=np.arange(-4, 4.1, 2),
        xminorlocator="null",
        yminorlocator="null",
        ylocator=(first.hlayers.values.astype(int)) / 100,
        suptitle=f"distribution of different periods with standardization of {std_type}",
    )

    modes = ["NAO", "EA"]

    for i, ax in enumerate(axes):
        df = xr2df(first, last, mode=modes[i])
        g = sns.violinplot(
            data=df,
            y="hlayers",
            x=modes[i],
            hue="period",
            kind="violin",
            palette="pastel",
            orient="h",
            ax=ax,
            split=split,
            dodge=True,
            linewidth=1,
        )
        g.axes.legend().remove()
        ax.set_xlim(-5, 5)
        ylabel = (first.hlayers.values) / 100
        ax.format(title=modes[i], xlabel="std", ylabel="gph/hpa")

    axes[-1].legend(loc="lr", ncols=1, title="period")
    return g
