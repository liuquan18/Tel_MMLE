#%%
import xarray as xr
import numpy as np
import pandas as pd
import random
import os
from scipy.signal import detrend

import src.Teleconnection.tools as tools
import src.Teleconnection.rolling_eof as rolling_eof
import src.Teleconnection.tools as tools
import warnings
import glob
# %%
from eofs.xarray import Eof
import src.Teleconnection.spatial_pattern as ssp

# %%
# read data
def read_EAR_JJA_data():
    odir = '/work/mh0033/m300883/Tel_MMLE/data/ERA5/'
    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + "zg_" + month + "/ts_era5_1940_2022.nc"

        # read data of single month
        zg_data = xr.open_dataset(zg_path)
        zg_data = zg_data.var129
        zg_data = zg_data/9.80665 # convert to geopotential height (unit m)

        # select 500hPa
        zg_data = zg_data.sel(plev = 50000)

        # detrend over every grid
        print("detrending the data ...")
        arr = zg_data.values
        arr = detrend(arr, axis=0, type='linear')
        zg_data = zg_data.copy(data=arr)

        data_JJA.append(zg_data)

    # combine JJA 
    data_JJA = xr.concat(data_JJA, dim="time")
    data_JJA = data_JJA.sortby("time")

    return data_JJA

#%%
def decompose_EAR_JJA(data):
    # using 
    data = data.transpose("time", "lat", "lon")
    wgts = tools.sqrtcoslat(data)
    # create solver
    solver = Eof(data, weights=wgts,center=True)
    # solve
    eof = solver.eofs(neofs=1)
    pc = solver.pcs(npcs=1)

    # deweight
    eof = eof / wgts[0]
    # 
    std_pc = pc.std(dim="time")
    pc = pc / std_pc
    std_pc = std_pc.squeeze()
    eof = eof * std_pc

    return eof,pc

# %%
data_JJA = read_EAR_JJA_data()
# %%
eof_result = ssp.doeof(data_JJA,nmode = 2,dim = 'time',standard = 'pc_temporal_std' )
# %%
eof_result = ssp.doeof(data_JJA,nmode = 2,dim = 'time',standard = 'eof_spatial_std' )

# %%
