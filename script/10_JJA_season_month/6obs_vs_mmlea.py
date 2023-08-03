# %%
import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import proplot as pplt
# %%
def read_eof_result(model):
    odir = '/work/mh0033/m300883/Tel_MMLE/data/ERA5/EOF_result/'
    filename = f'plev_50000_1940_2022_{model}_all.nc'
    ds = xr.open_dataset(odir+filename)
    return ds
# %%
