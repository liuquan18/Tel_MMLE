# %%
import xarray as xr
import numpy as np
import pandas as pd
import random
import os


import glob



# %%
def linear_trend(xarr):
    linear_coef = xarr.polyfit(dim="time", deg=1)
    linear_fitted = xr.polyval(xarr.time, linear_coef.polyfit_coefficients)
    return linear_fitted
def quadratic_trend(xarr):
    quadratic_coef = xarr.polyfit(dim="time", deg=2)
    quadratic_fitted = xr.polyval(xarr.time, quadratic_coef.polyfit_coefficients)
    return quadratic_fitted

def detrend(data,method = 'linear_trend'):
    ens_data = data.copy()
    try:
        ens_data = ens_data.mean(dim = 'ens')
    except ValueError:
        pass
    if method == 'linear_trend':
        fitted = ens_data.groupby('time.month').apply(linear_trend)
    elif method == 'quadratic_trend':
        fitted = ens_data.groupby('time.month').apply(quadratic_trend)
    detrended = data - fitted
    return detrended
