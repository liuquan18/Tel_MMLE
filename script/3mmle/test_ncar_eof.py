
#%%

import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools
import xarray as xr

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
#%%
# select trop
zg_ex = zg_demean.sel(hlayers = 20000)
# %%
zg_ex = zg_ex.stack(com = ('ens','time'))
# %%
eof,pc,fra = ssp.doeof(zg_ex)
# %%
