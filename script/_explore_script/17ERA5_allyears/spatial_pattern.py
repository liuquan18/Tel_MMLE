#%%
import xarray as xr
import numpy as np
import src.Teleconnection.spatial_pattern as ssp
from src.reanalysis.utils import read_gph_data

# %%
zg = read_gph_data("ERA5", plev = 50000, variable = 'zg', start_year = '1979',end_year = '2024')
# %%
zg = zg.load()
# %%
# %%
def decompose_period(xarr, nmode = 2,period = False):
    xarr = xarr.fillna(0) # fill nan with 0
    field = xarr.sortby("time")
    if 'ens' in xarr.dims:
        field = field.stack(com=("ens", "time"))
        dim = 'com'
    else: 
        dim = 'time'
    standard = 'eof_spatial_std' if period else 'pc_temporal_std'
    eof_result = ssp.doeof(field, standard=standard,nmode=nmode,dim = dim)
    return eof_result

# %%
eof_result = decompose_period(zg, nmode = 2, period = False)
# %%
eof_result.to_netcdf("/work/mh0033/m300883/High_frequecy_flow/data/ERA5/EOF_result/eof_result_Z500_1979_2024.nc")
# %%
