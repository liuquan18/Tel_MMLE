# %%
# import xarray, numpy 
import xarray as xr
import numpy as np
# %%
all_eof = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/ind_all_eof_result.nc")
# %%
all_pc = all_eof.pc
# %%
import matplotlib.pyplot as plt
ax = plt.subplot()
all_pc.sel(plev = 50000,mode = 'NAO').isel(time = slice(0,10)).plot.hist(ax = ax,color = 'b',alpha = 0.5) 
all_pc.sel(plev = 50000,mode = 'NAO').isel(time = slice(-10,None)).plot.hist(ax = ax,color = 'r',alpha = 0.5)
# %%
# standardize all_pc
all_pc_std = (all_pc - all_pc.mean(dim = 'time'))/all_pc.std(dim = 'time')


# %%
change_eof = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/ind_decade_eof_result.nc")
# %%
change_pc = change_eof.pc
# %%
