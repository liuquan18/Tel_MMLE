#%%
import importlib

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools

#%%
importlib.reload(ssp)

#%%
# genrate random data
time = pd.date_range("2010-10-01", "2020-10-01", freq="Y")
lon = np.linspace(-180, 180, 10)
lat = np.linspace(-90, 90, 10)
ens = [1, 2, 3]

np.random.seed(0)
values = np.random.random((10, 10, 10, 3))  # time, lon, lat, ens

data = xr.DataArray(
    values,
    dims=["time", "lon", "lat", "ens"],
    coords={"time": time, "lon": lon, "lat": lat, "ens": ens},
)


#%%
s_data = tools.stack_ens(data, withdim="time")

#%%
eof= ssp.doeof(s_data, nmode=2, dim="com")

#%%
# do eof
def test_doeof():
    assert eof.pc.std().values > 0 and eof.pc.std().values < 2


# %%
