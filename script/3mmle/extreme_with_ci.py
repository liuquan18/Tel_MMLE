# %%
# import xarray, numpy, pandas, proplot
import xarray as xr
import numpy as np
import pandas as pd
import proplot as pplt
import matplotlib.pyplot as plt

import src.extreme.extreme_ci as extreme

#%%
# Load data
ds = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/ind_first_pc.nc"
)
pc = ds.pc
# %%
firstpc = pc.sel(time=slice("1850", "1860"))
lastpc = pc.sel(time=slice("1990", "1999"))
# %%
firstcount = extreme.extreme_count_xr(firstpc, ci=True)
lastcount = extreme.extreme_count_xr(lastpc, ci=True)


# %%
