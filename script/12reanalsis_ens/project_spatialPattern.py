#%%
##### for projecting the spatial patterns using the index standardized by the first 10 years

##### SMILEs

import xarray as xr
import numpy as np
from scipy import stats
import pandas as pd
# %%
import src.MMLE_TEL.index_generator as index_generator
# %%
def read_eof_decade(model,fixed_pattern = 'decade_mpi'):
    """read eofs that is decomposed by decade"""
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_eof_result.nc"
    ds = xr.open_dataset(odir + filename)
    ds = ds.sel(mode="NAO")
    return ds

def read_zg_data(model):
    odir =  f"/work/mh0033/m300883/Tel_MMLE/data/{model}/"
    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + "zg_" + month + "/"
        data_JJA.append(index_generator.read_data(zg_path, plev = 50000))
    data = xr.concat(data_JJA, dim="time").sortby("time")
    return data

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
def project_period(zg_data, pc_data, period = 'first'):
    zg_data = zg_data.sortby('time')
    pc_data = pc_data.sortby('time')
    if period == 'first':
        p_field = zg_data.isel(time = slice(0,30)) # the first 10 years, 30 months
        p_pc = pc_data.sel(time = p_field.time, method = 'nearest')
    elif period == 'last':
        p_field = zg_data.sel(time = slice('2090','2099'))
        p_pc = pc_data.sel(time = slice('2090','2099'))
        if p_field.time.size == 0:
            p_field = zg_data.isel(time = slice(-30,None))
            p_pc = pc_data.sel(time = p_field.time, method = 'nearest')
    else:
        raise ValueError("period must be first or last")
    spatial_pattern = projected_pattern(p_field,p_pc)
    return spatial_pattern


# %%
def real_pattern(model):
    zg_data = read_zg_data(model)
    eof_data = read_eof_decade(model,fixed_pattern = 'decade_mpi')

    zg_data.load()
    eof_data.load()

    pc_data = eof_data.pc

    first_pattern = project_period(zg_data,pc_data,period = 'first')
    first_pattern.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/first_pattern_projected.nc")

    last_pattern = project_period(zg_data,pc_data,period = 'last')
    last_pattern.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/last_pattern_projected.nc")

# %%
import sys
models =  ["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"] #'MPI_GE"
mindex = int(sys.argv[1])

model = models[mindex - 1]
print(f"projecting the spatial pattern of {model} ...")
real_pattern(model)

