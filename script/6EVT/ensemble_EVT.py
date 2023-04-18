#%%
# imports
import numpy as np
import xarray as xr
import pandas as pd
from pyextremes import EVA
# %%
# load data
eof = xr.open_dataset('/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/ind_first_eof_result.nc')
# %%
pc = eof.pc
# %%
NAO500 = pc.sel(plev = 50000, mode = 'NAO')
# %%
# %%
ex = NAO500.copy()
ex = ex.isel(time = slice(0,10))
ex = ex.stack(com = ('time','ens'))
# %%
ex = ex.drop_vars(['time','ens']).squeeze()
# %%
# a new dim called time, which is a 1000-year-long time object
years = np.arange(2000,3001)
years = years.astype(str)
time = pd.to_datetime(years)
# %%
# generate a pd.date_range object, 1000 years long, with freq = 'Y'
