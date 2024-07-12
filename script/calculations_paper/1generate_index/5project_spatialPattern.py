#%%
##### for projecting the spatial patterns using the index standardized by the first 10 years

##### SMILEs

import xarray as xr
import numpy as np
from scipy import stats
import pandas as pd
from dateutil.relativedelta import relativedelta
import sys


# %%
import src.MMLE_TEL.index_generator as index_generator
# %%
def read_eof_plev_decade(model,fixed_pattern = 'decade_mpi'):
    """read eofs that is decomposed by decade"""
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_eof_result.nc"
    ds = xr.open_dataset(odir + filename)
    ds = ds.sel(mode="NAO")
    return ds

def read_eof_troposphere(model):
    eof = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_CMIP6/EOF_result/troposphere_ind_decade_first_JJA_eof_result.nc")
    eof = eof.sel(mode = 'NAO')
    return eof

def read_zg_data(model, plev =50000):
    odir =  f"/work/mh0033/m300883/Tel_MMLE/data/{model}/"
    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + "zg_" + month + "/"
        zg_month = index_generator.read_data(zg_path)
        if plev is not None:
            zg_month = zg_month.sel(plev = plev)
        data_JJA.append(zg_month)
    data = xr.concat(data_JJA, dim="time").sortby("time")
    return data

#%%
def _project(x,y):
    return stats.linregress(x,y)[0]

# %%

def projected_pattern(p_field, p_pc):
    p_field = p_field.stack(space = ['lat','lon'])
    p_field = p_field.stack(temporal = ['ens','time'])
    p_pc = p_pc.stack(temporal = ['ens','time'])
    p_field['temporal'] = np.arange(p_field.temporal.size)
    p_pc['temporal'] = np.arange(p_pc.temporal.size)

    spatial_pattern = xr.apply_ufunc(_project,
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
    try:
        zg_data = zg_data.sortby('plev')
        pc_data = pc_data.sortby('plev')
    except KeyError:
        print("plev is not found")
        pass
    if period == 'first':
        start_year = pc_data.time[0].dt.year.values
        end_year = start_year + 9

        p_field = zg_data.sel(time = slice(str(start_year), str(end_year))) # the first 10 years, 30 months
        p_pc = pc_data.sel(time = p_field.time, method = 'nearest')

    elif period == 'last':
        end_year = pc_data.time[-1].dt.year.values
        start_year = end_year - 9
        p_field = zg_data.sel(time = slice(str(start_year), str(end_year))) # the first 10 years, 30 months
        p_pc = pc_data.sel(time = p_field.time, method = 'nearest')

        if p_field.time.size == 0:
            p_field = zg_data.isel(time = slice(-30,None))
            p_pc = pc_data.sel(time = p_field.time, method = 'nearest')
    else:
        raise ValueError("period must be first or last")
    spatial_pattern = projected_pattern(p_field,p_pc)
    return spatial_pattern
#%%
def projet_all_period(zg_data,eof_data):
    intervals = eof_data.eof.decade

    real_patterns = []

    for start_year in intervals:
        start = pd.to_datetime(start_year.values)
        _end = start + relativedelta(years = 9)
        end = _end + relativedelta(months = 3)

        p_field = zg_data.sel(time = slice(start,end))
        p_pc = eof_data.pc.sel(time = slice(start,end))

        spatial_pattern = projected_pattern(p_field,p_pc)
        real_patterns.append(spatial_pattern)
    Real_patterns = xr.concat(real_patterns,dim = intervals)
    return Real_patterns

# %%
def real_pattern(zg_data, eof_data):

    zg_data.load()
    eof_data.load()

    pc_data = eof_data.pc

    first_pattern = project_period(zg_data,pc_data,period = 'first')
    first_pattern.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/first_pattern_projected.nc")

    last_pattern = project_period(zg_data,pc_data,period = 'last')
    last_pattern.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/last_pattern_projected.nc")

# %%
# for models in CMIP5
models =  ["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"] #'MPI_GE"
mindex = int(sys.argv[1])

model = models[mindex - 1]
print(f"projecting the spatial pattern of {model} ...")
zg_data = read_zg_data(model)
eof_data = read_eof_plev_decade(model,fixed_pattern = 'decade_mpi')

real_pattern(zg_data, eof_data)

#%%
# for MPI_GE_CMIP6
zg_data = read_zg_data("MPI_GE_CMIP6", plev = None)
eof_data = read_eof_troposphere("MPI_GE_CMIP6")
pc_data = eof_data.pc
# %%
zg_data.load()
pc_data.load()
#%%
zg_data = zg_data.sortby('plev')
pc_data = pc_data.sortby('plev')

# %%
first_pattern = project_period(zg_data,pc_data,period = 'first')

# %%
last_pattern = project_period(zg_data,pc_data,period = 'last')
# %%
first_pattern.to_netcdf("/work/mh0033/m300883/High_frequecy_flow/data/MPI_GE_CMIP6/season/EOF_result/troposphere_first_pattern_projected.nc")
# %%
last_pattern.to_netcdf("/work/mh0033/m300883/High_frequecy_flow/data/MPI_GE_CMIP6/season/EOF_result/troposphere_last_pattern_projected.nc")
# %%
