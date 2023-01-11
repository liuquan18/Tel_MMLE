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
mode = ["NAO","EA"]
values = np.array([[[np.arange(-5,6,1)]*10]*2]*2).reshape(11,10,2,2) # time,  ens, hlayers,mode

ex = xr.DataArray(
    values,
    dims=["time","ens", "hlayers","mode"],
    coords={"time": time, "ens": ens, "hlayers": hlayers,"mode":mode},
)
ex.name = 'pc'
# %%
pos, mpos, neg, mneg = EVT.mode_return_period(
    ex, mode = 'NAO',periods = [slice('2010','2015')],hlayers = 50000
)
# %%

def test_return_period():
    assert pos[0].shape[0] == 6