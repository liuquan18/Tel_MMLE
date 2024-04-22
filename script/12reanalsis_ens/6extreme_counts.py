# %%
import xarray as xr
import numpy as np
import pandas as pd
import os
# %%
import src.extreme.extreme_ci as extreme
#%%
import importlib
importlib.reload(extreme)
# %%
# read eof 
def read_eof(model, plev = 50000):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    first_eof_path = odir + f"EOF_result/first_plev{str(plev)}_eof.nc"
    last_eof_path = odir + f"EOF_result/last_plev{str(plev)}_eof.nc"

    first_eof = xr.open_dataset(first_eof_path)
    last_eof = xr.open_dataset(last_eof_path)
    return first_eof.pc, last_eof.pc

#%%
def extreme_count(model, plev = 50000):
    
    first_ind, last_ind = read_eof(model,plev = plev)
    
    first_extc = extreme.extreme_count_xr(first_ind,'bootstrap')
    last_extc = extreme.extreme_count_xr(last_ind,'bootstrap')

    # save
    save_path = f'/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/'
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    first_extc.to_netcdf(save_path + f'first_plev{str(plev)}_extc.nc')
    last_extc.to_netcdf(save_path + f'last_plev{str(plev)}_extc.nc')
    return first_extc, last_extc
# %%
plevs = [40000, 30000, 20000] #92500, 85000, 70000, 50000, 
for plev in plevs:
    print(f'plev{plev}')
    CR20_first_extc, CR20_last_extc = extreme_count('CR20_allens',plev)
# %%
# for "20CR"
    
plevs = [1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200]
model="CR20"

for plev in plevs:
    print(f'plev{plev}')
    CR20_first_extc, CR20_last_extc = extreme_count(model,plev)
# %%
