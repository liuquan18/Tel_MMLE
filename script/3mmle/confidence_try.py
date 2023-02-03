#%%
# import xarray, numpy, pandas, scipy
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.stats import bootstrap


# %%
# functions to do bootstrap on numpy array
def extreme_count(ts,threshold,gt=True):
    """
    count the number of extreme events in a time series
    ts: time series
    threshold: the threshold of extreme events
    gt: if True, count the events that are larger than the threshold
        if False, count the events that are smaller than the threshold
    """
    if gt:
        count = np.sum(ts>threshold)
    else:
        count = np.sum(ts<threshold)
    return count

def _pos_count(ts):
    return extreme_count(ts,2,True)
def _neg_count(ts):
    return extreme_count(ts,-2,False)

def bootstrap_pos_count_low(ts):
    res = bootstrap((ts,),_pos_count,random_state=0,vectorized=False)
    return res.confidence_interval.low
def bootstrap_pos_count_high(ts):
    res = bootstrap((ts,),_pos_count,random_state=0,vectorized=False)
    return res.confidence_interval.high

def bootstrap_neg_count_low(ts):
    res = bootstrap((ts,),_neg_count,random_state=0,vectorized=False)
    return res.confidence_interval.low
def bootstrap_neg_count_high(ts):
    res = bootstrap((ts,),_neg_count,random_state=0,vectorized=False)
    return res.confidence_interval.high


#%%
# apply the above functions on xarray DataArray
def extreme_count_xr(pc,ci = True):
    """
    calculate the extreme events of the stacked pc
    pc: xarray dataarray,with dimension ('hlayers','mode','ens','time)
    ci: if True, calculate the confidence interval of extreme events
    return: xarray dataarray
    """
    stacked_pc = pc.stack(stacked=['time','ens'])

    # positive extreme events
    pos_count_true = xr.apply_ufunc(_pos_count,stacked_pc,
                                   input_core_dims=[('stacked',)],
                                   output_core_dims=[[]],
                                   vectorize=True,
                                   dask='parallelized',
                                   exclude_dims=set(('stacked',)),
                                   output_dtypes=[int])
    # negative extreme events
    neg_count_true = xr.apply_ufunc(_neg_count,stacked_pc,
                                    input_core_dims=[('stacked',)],
                                    output_core_dims=[[]],
                                    vectorize=True,
                                    dask='parallelized',
                                    exclude_dims=set(('stacked',)),
                                    output_dtypes=[int])
    if ci:
        pos_count_low = xr.apply_ufunc(bootstrap_pos_count_low,stacked_pc,  
                                    input_core_dims=[('stacked',)],
                                    output_core_dims=[[]],
                                    vectorize=True,
                                    dask='parallelized',
                                    exclude_dims=set(('stacked',)),
                                    output_dtypes=[int])
        pos_count_high = xr.apply_ufunc(bootstrap_pos_count_high,stacked_pc,
                                        input_core_dims=[('stacked',)],
                                        output_core_dims=[[]],
                                        vectorize=True,
                                        dask='parallelized',
                                        exclude_dims=set(('stacked',)),
                                        output_dtypes=[int])
        pos_count = xr.concat([pos_count_low,pos_count_true,pos_count_high],dim='confidence')
        pos_count['confidence'] = ['low','true','high']


        neg_count_low = xr.apply_ufunc(bootstrap_neg_count_low,stacked_pc,
                                        input_core_dims=[('stacked',)],
                                        output_core_dims=[[]],
                                        vectorize=True,
                                        dask='parallelized',
                                        exclude_dims=set(('stacked',)),
                                        output_dtypes=[int])
        neg_count_high = xr.apply_ufunc(bootstrap_neg_count_high,stacked_pc,
                                        input_core_dims=[('stacked',)],
                                        output_core_dims=[[]],
                                        vectorize=True,
                                        dask='parallelized',
                                        exclude_dims=set(('stacked',)),
                                        output_dtypes=[int])
        neg_count = xr.concat([neg_count_low,neg_count_true,neg_count_high],dim='confidence')
        neg_count['confidence'] = ['low','true','high']

        extr_count = xr.concat([pos_count,neg_count],dim='extr_type')
        extr_count['extr_type'] = ['pos','neg']
        return extr_count
    else:
        extr_count = xr.concat([pos_count_true,neg_count_true],dim='extr_type')
        extr_count['extr_type'] = ['pos','neg']
        return extr_count

# %%
# test extreme_count_xr