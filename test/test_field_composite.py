#%%
import xarray as xr
import numpy as np
import src.composite.field_composite as fcomp
import pandas as pd

#%%
import importlib
importlib.reload(fcomp)

#%%
# genrate random temporal index
time = pd.date_range("2010-10-01", "2021-10-01", freq="Y")

values_index = np.arange(-5, 6, 1)
values_index = np.array([[values_index]*3]*2)

index = xr.DataArray(values_index, dims=["mode",'ens',"time"], coords={"mode":["NAO","EA"],"ens":[1,2,3],"time": time})

values_data = np.ones((11,10,10,3))
values_data[0] = values_data[0]*3
values_data[-1] = values_data[-1]*2
lon = np.arange(10)
lat = np.arange(10)
data = xr.DataArray(values_data,dims = ['time','lat','lon','ens'],coords = {'time':time,'lat':lat,'lon':lon,'ens':[1,2,3]})

# %%

tel_field_composite = fcomp.Tel_field_composite(index,data,threshold=4)
# %%
# test
def test_extreme():
    assert tel_field_composite.sel(extr_type = 'pos').mean().values == 2
    assert tel_field_composite.sel(extr_type = 'neg').mean().values == 3
# %%
