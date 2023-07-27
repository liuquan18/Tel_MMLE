# %%
import src.extreme.extreme_ci as extreme
import xarray as xr
import numpy as np
import importlib
import pandas as pd
#%%
importlib.reload(extreme)

#%%
# genrate random temporal index
time = pd.date_range("2010-10-01", "2021-10-01", freq="Y")
ens = np.arange(10)
plev = np.linspace(1000, 100, 10)
values = np.array([[np.arange(-1, 1.1, 0.2)] * 10] * 10).reshape(
    -1, 10, 10
)  # time,  ens, plev

ex = xr.DataArray(
    values,
    dims=["time", "ens", "plev"],
    coords={"time": time, "ens": ens, "plev": plev},
)

#%%
pc = ex
#%%
# test extreme.extreme_count_xr with pc
extr_count = extreme.extreme_count_xr(pc, ci='bootstrap')
# %%
extr_count = extreme.extreme_count_xr(pc, ci='AR1')

# %%
