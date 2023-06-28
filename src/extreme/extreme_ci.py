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
import statsmodels.api as sm
import random


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
def AR1_simulations_gen(ts):
    # Fit an AR1 model to the data
    model = sm.tsa.ARIMA(ts, order=(1, 0, 0)).fit()
    # Generate 9000 realizations of model_first, each 1000 long
    n_realizations = 5000
    n_obs = len(ts)
    simulations = np.empty((n_realizations, n_obs))
    for i in range(n_realizations):
        simulations[i, :] = model.simulate(nsimulations=n_obs)
    return simulations


# function to get the confidence interval of count of extreme events using AR1 model
def AR1_pos_count_low(ts):
    random.seed(1)
    AR1_simulations = AR1_simulations_gen(ts)
    counts = (AR1_simulations > 1.5).sum(axis=1)
    percentiles = np.percentile(counts, 5)
    return percentiles


# same as above, but for pso_count_high
def AR1_pos_count_high(ts):
    random.seed(2)
    AR1_simulations = AR1_simulations_gen(ts)
    counts = (AR1_simulations > 1.5).sum(axis=1)
    percentiles = np.percentile(counts, 95)
    return percentiles


# same as above, but for neg_count_low
def AR1_neg_count_low(ts):
    random.seed(3)
    AR1_simulations = AR1_simulations_gen(ts)
    counts = (AR1_simulations < -1.5).sum(axis=1)
    percentiles = np.percentile(counts, 5)
    return percentiles


# same as above, but for neg_count_high
def AR1_neg_count_high(ts):
    random.seed(4)
    AR1_simulations = AR1_simulations_gen(ts)
    counts = (AR1_simulations < -1.5).sum(axis=1)
    percentiles = np.percentile(counts, 95)
    return percentiles


#%%
# apply the above functions on xarray DataArray
def extreme_count_xr(pc, ci="AR1"):
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
    if ci == "bootstrap":
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

    # AR1 model
    if ci == "AR1":
        pos_count_low = xr.apply_ufunc(
            AR1_pos_count_low,
            stacked_pc,
            input_core_dims=[("stacked",)],
            output_core_dims=[[]],
            vectorize=True,
            dask="parallelized",
            exclude_dims=set(("stacked",)),
            output_dtypes=[int],
        )

        pos_count_high = xr.apply_ufunc(
            AR1_pos_count_high,
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
            AR1_neg_count_low,
            stacked_pc,
            input_core_dims=[("stacked",)],
            output_core_dims=[[]],
            vectorize=True,
            dask="parallelized",
            exclude_dims=set(("stacked",)),
            output_dtypes=[int],
        )

        neg_count_high = xr.apply_ufunc(
            AR1_neg_count_high,
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

# %%
def decadal_extrc_tsurf(index: xr.DataArray, ext_counts_xr = None, temp: xr.DataArray = None, plev=None,ci = 'AR1'):
    """
    extract the extreme count and the mean surface temperature every ten years.
    **Arguments**
        *index* the DataArray of pc index
        *temp* the -fldmean, -yearmean -ensmean surface temperauter.
    **Return**
        *ext_counts* the DataArray of extreme count, *time, *extr_type, *mode
        *t_surf_mean* the mean t_surface (the increase of the temperature)
        the time here use the first year of the decade.
    """
    print("counting the occureance of extremes and signifcance interval ...")
    if plev is not None:
        index = index.sel(plev=plev)



    # start time
    time_s = index.time[::10]
    # end time
    time_e = index.time[9::10]

    # create slice for each decade
    decade_slice = [slice(s, e) for s, e in zip(time_s, time_e)]

    ext_counts = []
    t_surf_mean = []


    for time in decade_slice:
        print(f" extreme counting in the decade of {time.start.dt.year.values} - {time.stop.dt.year.values}")

        period_pc = index.sel(time=time)
        # ensure that there are 10 years of data in period_pc
        if period_pc.time.size != 10:
            print(f" the length of the period is {len(period_pc.time)}, skip this period")
            # rasing a warning
            warnings.warn(f" the length of the period is {len(period_pc.time)}")
            break
        time_tag = period_pc.time[0] # for reference 

        # extreme count
        if ext_counts_xr is None:
            period_ext_count = extreme_count_xr(period_pc, ci=ci)
            period_ext_count['time'] = time_tag
            # set time as the new dimension
            period_ext_count = period_ext_count.expand_dims('time')
            ext_counts.append(period_ext_count)

        # tsurf
        if temp is not None:
            period_tm = temp.sel(time = period_pc.time, method = 'nearest')
            period_mean_t = period_tm.mean(dim="time")
            period_mean_t['time'] = time_tag
            period_mean_t = period_mean_t.expand_dims('time')
            t_surf_mean.append(period_mean_t)

    if ext_counts_xr is None:
        ext_counts_xr = xr.concat(ext_counts, dim='time')
    if temp is not None:
        t_surf_mean_xr = xr.concat(t_surf_mean, dim='time')
        
    return ext_counts_xr, t_surf_mean_xr

#%%
# plot function, x-axis is the extreme events count, y-axis is the pressure level
# vertical line, and fill_betweenx the confidence interval
def plot_extreme_count(ext_count, ax=None, label=None, colored=False):
    """
    plot the vertical profile of extreme event count
    """
    color = None
    style = None
    if colored:
        if label == "first10":
            color = "#1f77b4"
            style = color
        elif label == "last10":
            color = "#ff7f0e"
            style = color
    else:
        if label == "first10":
            color = "gray"
            style = "k-"
        elif label == "last10":
            color = "gray"
            style = "k--"

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
def extreme_count_profile(first_count, last_count, colored=False, **kwargs):
    """
    plot the extreme event count profile for the NAO and EA,
    and positive and negative extreme events
    """
    # parameters from kwargs
    xlim = kwargs.pop("xlim", None)
    
    fig = pplt.figure(
        # space=0,
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
        xlocator=20,
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
