#%%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import proplot as pplt
# %%
wb = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/wb_Aug/onepct_1850-1999_ens_0069.gph500_08.nc")
# %%
