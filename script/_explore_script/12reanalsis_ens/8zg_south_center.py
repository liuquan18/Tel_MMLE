#%%
from src.reanalysis.utils import read_gph_data as read_gph_data_reanalysis
from src.MMLE_TEL.index_generator import read_data as read_gph_data_MMLE

# %%
import xarray as xr
import numpy as np
# %%
# 500hpa geopotential height over the south center
CR20_zg = read_gph_data_reanalysis("CR20_allens",external_forcing=None,start_year = "1850")

#%%
MPI_GE_zg = []
for month in ["Jun", "Jul", "Aug"]:
    print(f"reading the gph data of {month} ...")
    zg_path = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_" + month + "/"
    data_month = read_gph_data_MMLE(zg_path, plev=50000,remove_ens_mean=False)
    MPI_GE_zg.append(data_month)
MPI_GE_zg = xr.concat(MPI_GE_zg, dim="time").sortby("time")
#%%
CR20_zg_south = CR20_zg.sel(lon = slice(-30,10), lat = slice(60,40)).mean(dim = ('lon','lat'))
# %%
MPI_GE_zg_south = MPI_GE_zg.sel(lon = slice(-30,10), lat = slice(60,40)).mean(dim = ('lon','lat'))
# %%
CR20_zg_south.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/ens_variability/CR20_zg_south_center.nc")

#%%
MPI_GE_zg_south.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/ens_variability/MPI_GE_zg_south_center.nc")
# %%
