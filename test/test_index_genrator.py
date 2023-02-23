#%%
# import
import xarray as xr
import numpy as np
import src.MMLE_TEL.index_generator as index_generate
import src.Teleconnection.spatial_pattern as ssp

#%%

# config
v_eof = 'ind' # vertical_eof
fpattern = 'first' # fixed pattern

#%%
# generate index
## CANESM2
canesm = index_generate.decompose_fixedPattern("CanESM2",v_eof,fpattern)
# %%


# the same data, but use the basic doeof function
data = canesm.data.sel(plev = 50000, time = slice('1951-01-01','1960-12-31'))

# %%
data = data.stack(com = ('time','ens'))

# %%
eof_result = ssp.doeof(data, nmode = 2)
# %%
# test
def test_index_generator():
    assert canesm.eof_result.sel(plev = 50000).fra.values == eof_result.fra.values
