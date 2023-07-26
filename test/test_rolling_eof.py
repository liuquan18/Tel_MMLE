#%%
import importlib

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import src.Teleconnection.rolling_eof as rolling_eof
import src.Teleconnection.tools as tools

#%%
importlib.reload(rolling_eof)

#%%
# genrate random data
time = pd.date_range("2000-10-01", "2100-10-01", freq="Y")
lon = np.linspace(-180, 180, 10)
lat = np.linspace(-90, 90, 10)
ens = np.arange(3)
plev = np.linspace(1000, 100, 2)
values = np.random.random((100, 10, 10, 3, 2))  # time, lon, lat, ens, plev

ex = xr.DataArray(
    values,
    dims=["time", "lon", "lat", "ens", "plev"],
    coords={"time": time, "lon": lon, "lat": lat, "ens": ens, "plev": plev},
)
ex = ex.sel(plev = 1000)
#%%
# rolling eof on generated dataset.
# eof= rolling_eof.rolling_eof(ex, nmode=2, fixed_pattern="all")

#%%
eof = rolling_eof.rolling_eof(ex, nmode=2,  fixed_pattern="decade")

#%%
# test
def test_rolling_eof():
    assert eof.pc.std().values > 0 and eof.pc.std().values < 2


# %%
