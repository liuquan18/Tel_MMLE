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


def xr2df(first10_all_period, last10_all_period, mode, compare):
    """
    prepare dataframe for seaborn.
    """
    first = to_dataframe(first10_all_period, mode=mode, name=mode)
    last = to_dataframe(last10_all_period, mode=mode, name=mode)
    if compare == "CO2":
        keys = ["first10", "last10"]
    elif compare == "temp":
        keys = ["0K", "4K"]
    periods = pd.concat([first, last], keys=keys, names=["period", "plev", "com"])
    periods = periods.reset_index().set_index("com")
    periods["plev"] = periods["plev"] / 100
    return periods


def plot_vilion(first, last, compare="CO2", split=False):
    fig = pplt.figure(space=0, refwidth="25em")
    axes = fig.subplots(nrows=1, ncols=2)
    axes.format(
        abc="a",
        abcloc="ul",
        xlocator=np.arange(-4, 4.1, 2),
        xminorlocator="null",
        yminorlocator="null",
        ylocator=(first.plev.values.astype(int)) / 100,
        suptitle=f"distribution of different periods",
    )

    modes = ["NAO", "EA"]

    for i, ax in enumerate(axes):
        df = xr2df(first, last, mode=modes[i], compare=compare)
        g = sns.violinplot(
            data=df,
            y="plev",
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
        ylabel = (first.plev.values) / 100
        ax.format(title=modes[i], xlabel="std", ylabel="gph/hpa")

    axes[-1].legend(loc="lr", ncols=1, title="period")
    return g
