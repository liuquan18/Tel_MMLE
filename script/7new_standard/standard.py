#%%
import xarray as xr
import numpy as np

# %%
decade_none = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/ind_decade_none_eof_result.nc")
# %%
#%% 
# stadardize decade_non with temporal mean and std
decade_std = (decade_none - decade_none.mean(dim = 'time'))/decade_none.std(dim = 'time')
# %%
# save to the same folder, replace the name with "ind_decade_own_eof_result"
decade_std.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/ind_decade_own_eof_result.nc")
# %%
