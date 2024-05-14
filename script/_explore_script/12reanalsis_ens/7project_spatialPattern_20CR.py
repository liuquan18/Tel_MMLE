#%%

import xarray as xr
import numpy as np
from scipy import stats
import pandas as pd
import glob
# %%
import src.MMLE_TEL.index_generator as index_generator

#%%
from src.reanalysis.utils import read_gph_data

#%%
def project(x,y):
    return stats.linregress(x,y)[0]

# %%

def projected_pattern(p_field, p_pc):
    p_field = p_field.stack(space = ['lat','lon'])
    p_field = p_field.stack(temporal = ['ens','time'])
    p_pc = p_pc.stack(temporal = ['ens','time'])
    p_field['temporal'] = np.arange(p_field.temporal.size)
    p_pc['temporal'] = np.arange(p_pc.temporal.size)

    spatial_pattern = xr.apply_ufunc(project,
                                    p_pc,
                                    p_field,
                                    input_core_dims=[['temporal'],['temporal']],
                                    vectorize = True)
    spatial_pattern = spatial_pattern.unstack('space')
    return spatial_pattern



# %%
first_eof = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_first_40_eof_std.nc")
last_eof = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_last_40_eof_std.nc")
# %%
first_pc = first_eof.pc.sel(mode = 'NAO')
last_pc = last_eof.pc.sel(mode = 'NAO') 
first_pc.load()
last_pc.load()

# %%
zg_data = read_gph_data('CR20_allens',start_year = '1850', end_year = '2015')

zg_data.load()

# %%
first_zg = zg_data.sel(time = slice('1850','1889'))
last_zg = zg_data.sel(time = slice('1976','2015'))
# %%
first_pattern = projected_pattern(first_zg, first_pc)
first_pattern.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_first_40_pattern_projected.nc") 

#%%
last_pattern = projected_pattern(last_zg, last_pc)
last_pattern.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_last_40_pattern_projected.nc")
# %%
