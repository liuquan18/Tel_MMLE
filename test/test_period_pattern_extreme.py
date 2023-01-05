#%%
import xarray as xr
import pandas as pd
import numpy as np

import src.extreme.period_pattern_extreme as period_ext

#%%
import importlib

importlib.reload(period_ext)


#%%
# genrate random temporal index
period_time = pd.date_range("2010-10-01", "2021-10-01", freq="Y")
all_time = pd.date_range("2010-10-01", "2060-10-01", freq="Y")

period_values = np.arange(-5, 6, 1)
np.random.seed(0)
all_values = np.random.randn(len(all_time))

period_period = xr.DataArray(period_values, dims=["time"], coords={"time": period_time})

all_all = xr.DataArray(all_values, dims=["time"], coords={"time": all_time})


# %%
period_extreme_count = period_ext.period_extreme_count(
    period_period, all_all, dim="time", threshold=3
)
# %%
# test
def test_period_pattern_extreme():
    assert period_extreme_count.values.mean() == 2

# %%
