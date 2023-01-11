#%%
import src.EVT.return_period as EVT
import pandas as pd
import numpy as np
import xarray as xr

#%%
# genrate random temporal index
time = pd.date_range("2010-10-01", "2021-10-01", freq="Y")
ens = np.arange(10)
hlayers = [100000,50000]
values = np.array([[np.arange(-5,6,1)]*10]*2).reshape(11,10,2) # time,  ens, hlayers

ex = xr.DataArray(
    values,
    dims=["time","ens", "hlayers"],
    coords={"time": time, "ens": ens, "hlayers": hlayers},
)
# %%
