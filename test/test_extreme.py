#%%
import xarray as xr
import numpy as np
import pandas as pd
import pytest
# %%

import src.extreme.extreme as ext

#%%
import importlib
importlib.reload(ext)

#%%
# genrate random temporal index
time = pd.date_range("2010-10-01", "2021-10-01", freq="Y")
ens = np.arange(10)
hlayers = np.linspace(1000, 100, 10)
values = np.array([[np.arange(-5,6,1)]*10]*10).reshape(11,10,10) # time,  ens, hlayers

ex = xr.DataArray(
    values,
    dims=["time","ens", "hlayers"],
    coords={"time": time, "ens": ens, "hlayers": hlayers},
)
# %%
extreme = ext.extreme(ex,threshod = 4)
extreme_count = ext.count_extreme(extreme,dim = ('time','ens'))
# %%
# test
def test_extreme():
    assert extreme_count.mean().values == 10

# %%
