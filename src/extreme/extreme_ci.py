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
import warnings

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
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try: # in case of count=0, the bootstrap will give a warning
            res = bootstrap(
                (ts,), _pos_count, random_state=0, vectorized=False, confidence_level=cl
            )
            low = res.confidence_interval.low
        except Warning:
            low = 0
    return low


def bootstrap_pos_count_high(ts, cl=0.95):
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            res = bootstrap(
                (ts,), _pos_count, random_state=0, vectorized=False, confidence_level=cl
            )
            high = res.confidence_interval.high
        except Warning:
            high = 0
    return high


def bootstrap_neg_count_low(ts, cl=0.95):
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            res = bootstrap(
                (ts,), _neg_count, random_state=0, vectorized=False, confidence_level=cl
            )
            low = res.confidence_interval.low
        except Warning:
            low = 0
    return low



def bootstrap_neg_count_high(ts, cl=0.95):
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            res = bootstrap(
                (ts,), _neg_count, random_state=0, vectorized=False, confidence_level=cl
            )
            high = res.confidence_interval.high
        except Warning:
            high = 0
    return high



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
def decadal_extrc(index: xr.DataArray, plev=None,ci = 'bootstrap',window = 10):
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
    years = np.unique(index.time.dt.year).astype('str')
    time_s = years[::window]
    # end time
    time_e = years[window-1::window]

    # create slice for each decade
    decade_slice = [slice(s, e) for s, e in zip(time_s, time_e)]

    ext_counts = []
    for time in decade_slice:
        print(f" extreme counting in the decade of {time.start} - {time.stop}")

        period_pc = index.sel(time=time)
        # ensure that there are 10 years of data in period_pc
        if len(np.unique(period_pc.time.dt.year)) != 10:
            print(f" the length of the period is {len(period_pc.time)}, skip this period")
            # rasing a warning
            warnings.warn(f" the length of the period is {len(period_pc.time)}")
            break
        time_tag = period_pc.time[0] # for reference 

        # extreme count
        period_ext_count = extreme_count_xr(period_pc, ci=ci)
        period_ext_count['time'] = time_tag
        # set time as the new dimension
        period_ext_count = period_ext_count.expand_dims('time')
        ext_counts.append(period_ext_count)
        ext_counts_xr = xr.concat(ext_counts, dim='time')
        
    return ext_counts_xr

#%%
def decade_tsurf(tsurf):

    start = tsurf.time.dt.year.values[0]
    end = tsurf.time.dt.year.values[-1]
    
    time_s = np.arange(start,end,10)
    time_e = np.arange(start+9,end+1,10)
    
    time_tag = pd.date_range(str(time_s[0]) + '-06-16',periods=len(time_s),freq = '10Y')
    
    decade_slices = [slice(str(s), str(e)) for s, e in zip(time_s, time_e)]

    tsurf_dec_mean = [
    tsurf.sel(time=decade_slice).mean(dim="time") for decade_slice in decade_slices
    ]
    tsurf_dec_mean = xr.concat(tsurf_dec_mean, dim = time_tag).rename({'concat_dim':'time'})
    return tsurf_dec_mean.squeeze()
