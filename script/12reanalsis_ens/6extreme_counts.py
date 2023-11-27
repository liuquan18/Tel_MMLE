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
def read_eof(model,group_size = 40):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    first_eof_path = odir + f"EOF_result/first_{str(group_size)}_eof_std.nc"
    last_eof_path = odir + f"EOF_result/last_{str(group_size)}_eof_std.nc"

    first_eof = xr.open_dataset(first_eof_path)
    last_eof = xr.open_dataset(last_eof_path)
    return first_eof.pc, last_eof.pc

#%%
def extreme_count(model,group_size):
    
    first_ind, last_ind = read_eof(model,group_size=group_size)
    
    first_extc = extreme.extreme_count_xr(first_ind,'bootstrap')
    last_extc = extreme.extreme_count_xr(last_ind,'bootstrap')

    # save
    save_path = f'/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/'
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    first_extc.to_netcdf(save_path + f'first_{str(group_size)}_extc.nc')
    last_extc.to_netcdf(save_path + f'last_{str(group_size)}_extc.nc')
    return first_extc, last_extc
# %%
CR20_first_extc, CR20_last_extc = extreme_count('CR20',40)

# %%
