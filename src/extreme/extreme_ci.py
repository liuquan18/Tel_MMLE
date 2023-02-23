#%%
# import xarray, numpy, pandas, scipy
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import scipy.stats as stats
import matplotlib.pyplot as plt
import proplot as pplt
from scipy.stats import bootstrap


# %%
# functions to do bootstrap on numpy array
def extreme_count(ts, threshold, gt=True):
    """
    count the number of extreme events in a time series
    ts: time series
    threshold: the threshold of extreme events
    gt: if True, count the events that are larger than the threshold
        if False, count the events that are smaller than the threshold
    """
    if gt:
        count = np.sum(ts > threshold)
    else:
        count = np.sum(ts < threshold)
    return count


def _pos_count(ts):
    return extreme_count(ts, 1.5, True)


def _neg_count(ts):
    return extreme_count(ts, -1.5, False)


def bootstrap_pos_count_low(ts, cl=0.95):
    res = bootstrap(
        (ts,), _pos_count, random_state=0, vectorized=False, confidence_level=cl
    )
    return res.confidence_interval.low


def bootstrap_pos_count_high(ts, cl=0.95):
    res = bootstrap(
        (ts,), _pos_count, random_state=0, vectorized=False, confidence_level=cl
    )
    return res.confidence_interval.high


def bootstrap_neg_count_low(ts, cl=0.95):
    res = bootstrap(
        (ts,), _neg_count, random_state=0, vectorized=False, confidence_level=cl
    )
    return res.confidence_interval.low


def bootstrap_neg_count_high(ts, cl=0.95):
    res = bootstrap(
        (ts,), _neg_count, random_state=0, vectorized=False, confidence_level=cl
    )
    return res.confidence_interval.high


#%%
# apply the above functions on xarray DataArray
def extreme_count_xr(pc, ci=True):
    """
    calculate the extreme events of the stacked pc
    pc: xarray dataarray,with dimension ('plev','mode','ens','time)
    ci: if True, calculate the confidence interval of extreme events
    return: xarray dataarray
    """
    stacked_pc = pc.stack(stacked=["time", "ens"])

    # positive extreme events
    pos_count_true = xr.apply_ufunc(
        _pos_count,
        stacked_pc,
        input_core_dims=[("stacked",)],
        output_core_dims=[[]],
        vectorize=True,
        dask="parallelized",
        exclude_dims=set(("stacked",)),
        output_dtypes=[int],
    )
    # negative extreme events
    neg_count_true = xr.apply_ufunc(
        _neg_count,
        stacked_pc,
        input_core_dims=[("stacked",)],
        output_core_dims=[[]],
        vectorize=True,
        dask="parallelized",
        exclude_dims=set(("stacked",)),
        output_dtypes=[int],
    )
    if ci:
        pos_count_low = xr.apply_ufunc(
            bootstrap_pos_count_low,
            stacked_pc,
            input_core_dims=[("stacked",)],
            output_core_dims=[[]],
            vectorize=True,
            dask="parallelized",
            exclude_dims=set(("stacked",)),
            output_dtypes=[int],
        )
        pos_count_high = xr.apply_ufunc(
            bootstrap_pos_count_high,
            stacked_pc,
            input_core_dims=[("stacked",)],
            output_core_dims=[[]],
            vectorize=True,
            dask="parallelized",
            exclude_dims=set(("stacked",)),
            output_dtypes=[int],
        )
        pos_count = xr.concat(
            [pos_count_low, pos_count_true, pos_count_high], dim="confidence"
        )
        pos_count["confidence"] = ["low", "true", "high"]

        neg_count_low = xr.apply_ufunc(
            bootstrap_neg_count_low,
            stacked_pc,
            input_core_dims=[("stacked",)],
            output_core_dims=[[]],
            vectorize=True,
            dask="parallelized",
            exclude_dims=set(("stacked",)),
            output_dtypes=[int],
        )
        neg_count_high = xr.apply_ufunc(
            bootstrap_neg_count_high,
            stacked_pc,
            input_core_dims=[("stacked",)],
            output_core_dims=[[]],
            vectorize=True,
            dask="parallelized",
            exclude_dims=set(("stacked",)),
            output_dtypes=[int],
        )
        neg_count = xr.concat(
            [neg_count_low, neg_count_true, neg_count_high], dim="confidence"
        )
        neg_count["confidence"] = ["low", "true", "high"]

        extr_count = xr.concat([pos_count, neg_count], dim="extr_type")
        extr_count["extr_type"] = ["pos", "neg"]
        return extr_count
    else:
        extr_count = xr.concat([pos_count_true, neg_count_true], dim="extr_type")
        extr_count["extr_type"] = ["pos", "neg"]
        return extr_count


#%%
# plot function, x-axis is the extreme events count, y-axis is the pressure level
# vertical line, and fill_betweenx the confidence interval
def plot_extreme_count(ext_count, ax=None, label=None, colored=False):
    """
    plot the vertical profile of extreme event count
    """
    if colored:
        if label == "first10":
            style = color = "#1f77b4"
        elif label == "last10":
            style = color = "#ff7f0e"
    else:
        if label == "first10":
            style = "k-"
            color = "gray"
        elif label == "last10":
            style = "k--"
            color = "gray"

    if ax is None:
        ax = plt.gca()

    y = ext_count.plev / 100
    true = ext_count.sel(confidence="true").values
    low = ext_count.sel(confidence="low").values
    high = ext_count.sel(confidence="high").values

    # plot the confidence interval
    ax.fill_betweenx(y, low, high, color=color, alpha=0.3)
    # plot the true value
    line = ax.plot(true, y, style, linewidth=1, label=label)
    return line


# %%
def extreme_count_profile(
    first_count,
    last_count,
    colored=False,
    xlim = (-5, 45)
):
    """
    plot the extreme event count profile for the NAO and EA,
    and positive and negative extreme events
    """
    fig = pplt.figure(
        space=0,
        refwidth="20em",
    )
    axes = fig.subplots(nrows=2, ncols=2)
    axes.format(
        abc=True,
        abcloc="ul",
        abcstyle="a",
        xlabel="Extreme event count",
        ylabel="Pressure level (hPa)",
        suptitle="Extreme event count profile",
        ylim=(1000, 200),
        xminorticks="null",
        yminorticks="null",
        grid=False,
        toplabels=("pos", "neg"),
        leftlabels=("NAO", "EA"),
        xlocator=10,
    )

    # plot the extreme event count profile
    labels = ["first10", "last10"]
    # the default color of matplotlib

    for i, extreme_count in enumerate([first_count, last_count]):
        plot_extreme_count(
            extreme_count.sel(mode="NAO", extr_type="pos"),
            axes[0, 0],
            label=labels[i],
            colored=colored,
        )
        plot_extreme_count(
            extreme_count.sel(mode="NAO", extr_type="neg"),
            axes[0, 1],
            label=labels[i],
            colored=colored,
        )

        plot_extreme_count(
            extreme_count.sel(mode="EA", extr_type="pos"),
            axes[1, 0],
            label=labels[i],
            colored=colored,
        )
        plot_extreme_count(
            extreme_count.sel(mode="EA", extr_type="neg"),
            axes[1, 1],
            label=labels[i],
            colored=colored,
        )
    for ax in axes:
        ax.set_xlim(xlim)
    # add legend
    axes[0, 0].legend(loc="lr", ncols=1, frame=True)
