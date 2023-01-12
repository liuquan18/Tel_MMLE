#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.MMLE_TEL.index_generator as index_generate
import src.MMLE_TEL.quick_plot as quick_plot

# config
v_eof = 'ind' # vertical_eof
fpattern = 'first' # fixed pattern

# %%
canesm2 = index_generate.decompose_fixedPattern("test_CanESM2",v_eof,fpattern)

# %%
cesm = index_generate.decompose_fixedPattern("test_CESM",v_eof,fpattern,'temp')
# %%



#%%

import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools

#%%
import importlib

importlib.reload(ssp)

# %%
zg_path = "/work/mh0033/m300883/Tel_MMLE/data/test_CESM/zg_processed/*.nc"
zg_data = xr.open_mfdataset(zg_path, combine="nested", concat_dim="ens")


# %%
zg_data = zg_data.zg
zg_data = zg_data.rename({"plev": "hlayers"})  # to adapt

zg_ens_mean = zg_data.mean(dim="ens")
zg_demean = zg_data - zg_ens_mean

# select trop
zg_trop = zg_demean.sel(hlayers=slice(100000, 20000))
# %%
