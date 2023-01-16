#%%
import importlib

import numpy as np
import pandas as pd
import pytest
import xarray as xr

import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools
import src.MMLE_TEL.index_generator as index_generate

# %%
data = xr.open_mfdataset("/work/mh0033/m300883/Tel_MMLE/data/CESM1_CAM5/zg_processed/*.nc",concat_dim = 'ens', combine = 'nested')
# %%
ddata = data.isel(time = slice(0,10))
# %%
ddata = ddata.sel(plev = slice(100000,20000)).zg
# %%
cdata = ddata.stack(com = ("time","ens"))
# %%
eof,pc,fra = ssp.doeof(cdata.sel(plev = 50000))


# %%
cesm = index_generate.decompose_fixedPattern("test_CESM1_CAM5",'ind','first')

# %%



