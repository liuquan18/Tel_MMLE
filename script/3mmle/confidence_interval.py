#%% 
# import xarray, numpy, pandas and matplotlib
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#%%
dir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/ind_first_pc.nc"
ds = xr.open_dataset(dir)

# %%
pc = ds.pc
# %%
