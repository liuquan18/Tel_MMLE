# %%
import xarray as xr
import numpy as np
import pandas as pd
import random
import os

import src.Teleconnection.vertical_eof as vertical_eof
import src.Teleconnection.tools as tools
import src.Teleconnection.rolling_eof as rolling_eof
import src.warming_stage.warming_stage as warming_stage
import src.Teleconnection.spatial_pattern as ssp

import warnings
import glob

import proplot as pplt
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview

#%%
model = 'ERA5_allens'
fixedPattern="decade"
plev=50000
standard="first"
#%%
# read gph data
odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
ts_mean_path = odir + "ts_processed/ens_fld_year_mean.nc"
save_path = odir + "EOF_result/"

data_JJA = []
for month in ["Jun", "Jul", "Aug"]:
    print(f"reading the gph data of {month} ...")
    zg_path = odir + "zg_" + month + "/"
    data_JJA.append(read_data(zg_path, plev=plev, remove_ens_mean=False))
data = xr.concat(data_JJA, dim="time").sortby("time")

#%%
ens_mean = data.mean(dim = 'ens')
#%%
ens_mean_c = ens_mean.copy()
ens_mean_c['time'] = np.arange(ens_mean_c.time.size)
#%%
linear_coef = ens_mean_c.polyfit(dim='time', deg=1)
#%%
linear_fitted = xr.polyval(ens_mean_c.time, linear_coef.polyfit_coefficients)

linear_fitted['time'] = ens_mean.time
#%%
linear_detrend = data - linear_fitted

#%%
first_40 = linear_detrend.sel(time=slice('1940','1980'))

last_40 = linear_detrend.sel(time=slice('1982','2022'))

#%%
def decompose_40(xarr):
    field = xarr.sortby('time')
    field = field.stack(com = ('ens','time'))
    eof_result = ssp.doeof(field)
    return eof_result

#%%
first_40_eof = decompose_40(first_40)
last_40_eof = decompose_40(last_40)

#%%

first_40_eof['mode'] = ['None','NAO']
last_40_eof['mode'] = ['None','NAO']
#%%
def standard_40(first_40_eof,last_40_eof):
    mean = first_40_eof['pc'].mean(dim=('time','ens'))
    std = first_40_eof['pc'].std(dim=('time','ens'))

    first_40_eof['pc'] = (first_40_eof['pc'] - mean) / std
    last_40_eof['pc'] = (last_40_eof['pc'] - mean) / std
    return first_40_eof, last_40_eof

#%%

first_40_eof_std, last_40_eof_std = standard_40(first_40_eof,last_40_eof)

#%%
first_40_eof.to_netcdf(save_path + 'first_40_eof.nc')
last_40_eof.to_netcdf(save_path + 'last_40_eof.nc')

first_40_eof_std.to_netcdf(save_path + 'first_40_eof_std.nc')
last_40_eof_std.to_netcdf(save_path + 'last_40_eof_std.nc')
# %%
fig2 = pplt.figure(figsize=(180 / 25.4, 90 / 25.4),sharex=False,sharey=False)
fig2.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
    
)
gs = pplt.GridSpec(
    ncols=2,
    nrows=1,
)
ax1 = fig2.add_subplot(gs[0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax2 = fig2.add_subplot(gs[1])


ax1, fmap, lmap = stat_overview.spatial_pattern_plot(
    ax1,
    first_40_eof.eof.sel(mode = 'NAO',decade = '1940').squeeze(),
    first_40_eof.fra.sel(mode = 'NAO',decade = '1940').squeeze(),
    last_40_eof.eof.sel(mode = 'NAO',decade = '1982').squeeze(),
    last_40_eof.fra.sel(mode = 'NAO',decade = '1982').squeeze(),
    levels = np.arange(-2,2.1,0.4),
)
# %%
