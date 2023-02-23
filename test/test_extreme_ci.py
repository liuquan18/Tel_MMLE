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
values = np.array([[np.arange(-5, 6, 1)] * 10] * 10).reshape(
    11, 10, 10
)  # time,  ens, plev

ex = xr.DataArray(
    values,
    dims=["time", "ens", "plev"],
    coords={"time": time, "ens": ens, "plev": plev},
)

#%%
pc = ex
# test extreme.extreme_count_xr with pc
extr_count = extreme.extreme_count_xr(pc, ci=False, threshold=2)
# %%
pos_low = extreme.bootstrap_pos_count_low(pc, threshold=2)
# %%
