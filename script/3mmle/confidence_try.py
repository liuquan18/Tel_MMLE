#%%
# import xarray, numpy, pandas, scipy
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import scipy.stats as stats
import matplotlib.pyplot as plt
import src.extreme.period_pattern_extreme as extreme

import importlib
importlib.reload(extreme)


#%%
# create random xarray dataarray with dims of (time:50, ensemble:10),time starting from 1850-12-31
# and ensemble starting from 1
# dataarray is a random normal distribution with mean 0 and std 1
# dataarray is a 2d array with 50 rows and 10 columns
ds = xr.DataArray(np.random.randn(10,5), dims=('time','ens'), coords={'time':pd.date_range('1850-12-31', periods=10, freq='M'),'ens':np.arange(1,6)})


# %%
extr = extreme.period_extreme(ds,2)
# %%
count= extreme.extreme_count(extr, dim = ('time','ens'))
# %%
