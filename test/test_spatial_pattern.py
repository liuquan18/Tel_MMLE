#%%
import importlib

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools

#%%
import importlib

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
eof, pc, fra = ssp.doeof(s_data, nmode=2, dim="com", standard=True)

#%%
ppc = ssp.project_field(s_data, eof, dim="com", standard=True)

#%%
# do eof
def test_doeof():
    assert pc.std().values > 0 and pc.std().values < 2
    assert ppc.std().values > 0 and ppc.std().values < 2


# %%
