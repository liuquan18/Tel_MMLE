#%%
# imports
import numpy as np
import xarray as xr
import pandas as pd
from pyextremes import EVA
from datetime import date
import matplotlib.pyplot as plt
# %%
# load data
eof = xr.open_dataset('/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/ind_first_eof_result.nc')
# %%
pc = eof.pc
# %%
NAO500 = pc.sel(plev = 50000, mode = 'NAO')

# %%
def EVT_return_period(df,axes,extreme_type = 'hight',threshold = 1.5):
    model = EVA(df)

    model.get_extremes(method='POT', threshold=threshold,r = '1H',extremes_type=extreme_type)
    model.fit_model()
    model.plot_return_values(
                         alpha = 0.95,return_period_size='1H',ax=axes)
    # plt.ylim(1.5,3)
    # plt.xlim(10,100)
    # don't set log scale for x-axis
    plt.xscale('linear')
    return model

def post_process(ex):
    ex = ex.stack(com = ('time','ens'))
    ex = ex.drop_vars(['time','ens']).squeeze()
    # generate 1000 hours long pd.time_range object
    hours = pd.date_range(start = date(2000,1,1), periods = 1000, freq = 'H')
    # generate a pd.date_range object, 1000 years long, with freq = 'Y'
    df = pd.Series(ex.values, index = hours)
    return df
# %%
ex = NAO500.copy()
first = ex.isel(time = slice(0, 10))
first = post_process(first)

last = ex.isel(time = slice(-10, None))
last = post_process(last)

#%%
axes = plt.subplot()
# model = EVT_return_period(first,axes,extreme_type = 'high',threshold = 2.2)
model = EVT_return_period(last,axes,extreme_type = 'high',threshold = 2.2)
# plt.xlim(0,400)

#%%
axes = plt.subplot()
model = EVT_return_period(first,axes,extreme_type = 'low',threshold = -2)
model = EVT_return_period(last,axes,extreme_type = 'low',threshold = -2)



#%%
model = EVA(first)

model.get_extremes(method='POT', threshold=-2,r = '1H',extremes_type='low')
# model.fit_model()

#%%
model.plot_return_values(
                        alpha = None,return_period_size='1H',ax=axes)
# plt.ylim(1.5,3)
# plt.xlim(10,100)
# don't set log scale for x-axis
plt.xscale('linear')