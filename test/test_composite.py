#%%
import xarray as xr
import numpy as np
import pandas as pd
import src.composite.composite as comp
import pytest

# %%
import importlib
importlib.reload(comp)
#%%
# genrate random temporal index
time = pd.date_range("2010-10-01", "2021-10-01", freq="Y")

# random value between -3, 3
np.random.seed(0)
values_index = np.random.rand(11, 3) * 6 - 3

#%%
index = xr.DataArray(values_index, dims=["time","ens"], coords={"time": time, "ens": range(3)})

#%%
values_data = np.ones((11,10,10))
values_data[0] = values_data[0]*2
values_data[-1] = values_data[-1]*2
values_data = [values_data]*3
values_data = np.array(values_data).reshape(11,3,10,10)
lon = np.arange(10)
lat = np.arange(10)
data = xr.DataArray(values_data,dims = ['time','ens','lat','lon'],coords = {'time':time,'ens':np.arange(3),'lat':lat,'lon':lon})

#%%
index_c = index.copy()
data_c = data.copy()
index_c = index_c.stack(com = ['time','ens'])
data_c = data_c.stack(com = ['time','ens'])

#%%
composite_mean = comp.reduce_var(index_c,data_c,dim = 'com',reduction = 'mean',bootstrap = True)

#%%

first_index = index.isel(time = slice(0,5))
last_index = index.isel(time = slice(6,11))
#%%
first_index = first_index.expand_dims({'mode':['NAO']})
last_index = last_index.expand_dims({'mode':['NAO']})

# %%
composite_mean = comp.first_last_extreme_composite(
    first_index, last_index, data, dim="com", reduction="mean_same_number", bootstrap=True,count = 2
)
# %%