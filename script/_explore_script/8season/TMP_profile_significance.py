#%%
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
# %%
eof = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/troposphere_ind_decade_first_JJAS_eof_result.nc")
index = eof.pc
# %%
first = index.isel(time = slice(0, 10))
last = index.isel(time = slice(-10, None))
# %%
def bootstrap_diff_mean(mode, plev, extr_type, n_bootstrap=10000):
    # Define the significance level and number of bootstrap iterations
    alpha = 0.05
    n_bootstrap = n_bootstrap

    plevs = first.plev.values
    sigs = np.zeros(len(plevs))



#%%
def sig_test(first, last, mode, extr_type = 'pos',n_bootstrap = 10000, alpha = 0.05):

    # for all plevs
    plevs = first.plev.values
    sigs = np.zeros(len(plevs))

    for i, plev in enumerate(plevs):
        first_df = to_dataframe(first,mode, plev )
        last_df = to_dataframe(last,mode, plev )
        sigs[i] = boot(first_df, last_df, n_bootstrap, alpha, extr_type)

    return sigs

#%%
def boot(first_df, last_df, n_bootstrap = 10000, alpha = 0.05, extr_type = 'pos'):

    bootstrap = np.zeros(n_bootstrap)
    n = len(first_df)
    for i in range(n_bootstrap):
        sampy1 = first_df.sample(n, replace=True)
        sampy2 = last_df.sample(n, replace=True)

        count_first = extreme_count(sampy1.pc, 1.5, extr_type == 'pos')
        count_last = extreme_count(sampy2.pc, 1.5, extr_type == 'pos')
        
        bootstrap[i] = count_first - count_last
    # Calculate the confidence interval
    confidence = np.percentile(bootstrap, [100*alpha/2, 100-alpha/2])

    # Determine if the confidence interval covers zero
    significance = (confidence[0] > 0) or (confidence[1] < 0)
    return significance


#%%
def to_dataframe(xarr,mode,plev):
    xarr = xarr.sel(mode = mode, plev = plev).stack(com = ('ens', 'time'))
    xarr = xarr.drop(('ens', 'time'))
    df = xarr.to_dataframe()
    df = df.reset_index()
    df = df.drop(columns = ['mode', 'plev'])
    return df[['pc']]

# %%
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

# %%
