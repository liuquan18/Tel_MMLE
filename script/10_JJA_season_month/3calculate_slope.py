#%%
import xarray as xr
import numpy as np
import pandas as pd
import statsmodels.api as sm
import os

#%% calculate the slope of the extreme counts
def cal_slope(extrc):
    x = np.arange(len(extrc))
    x = sm.add_constant(x)
    y = extrc
    model = sm.OLS(y,x)
    results = model.fit()
    slope = results.params[1]
    return slope

def cal_conf_low(extrc):
    # calculate the 95% confidence interval
    x = np.arange(len(extrc))
    x = sm.add_constant(x)
    y = extrc
    model = sm.OLS(y,x)
    results = model.fit()
    confint = results.conf_int(alpha = 0.05, cols = None)[1][0] # the lower bound
    return confint

def cal_conf_high(extrc):
    x = np.arange(len(extrc))
    x = sm.add_constant(x)
    y = extrc
    model = sm.OLS(y,x)
    results = model.fit()
    confint = results.conf_int(alpha = 0.05, cols = None)[1][1] # the upper bound
    return confint

#%%
def read_extrc(model, fixed_pattern="decade_mpi"):
    """read extreme counts"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_count/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_extre_counts.nc"
    ds = xr.open_dataset(odir + filename).pc

    # divide the ensemble size of each model
    ens_sizes = {'MPI_GE':100,
             'MPI_GE_onepct':100,
             'CanESM2':50,
                'CESM1_CAM5':40,
                'MK36':30,
                'GFDL_CM3':20,
             }
    ds = ds/ens_sizes[model]
    return ds
#%%
def SMILE_slope(extrc,fixed_pattern = 'decade_mpi'):
    try:
        extrc = extrc.sel(confidence = 'true').drop('confidence')
    except KeyError:
        pass
    slope = xr.apply_ufunc(cal_slope,extrc,input_core_dims=[['time']],vectorize=True, dask='parallelized')
    conf_low = xr.apply_ufunc(cal_conf_low,extrc,input_core_dims=[['time']],vectorize=True, dask='parallelized')
    conf_high = xr.apply_ufunc(cal_conf_high,extrc,input_core_dims=[['time']],vectorize=True, dask='parallelized')
    new_dim = xr.IndexVariable("slopes", ['true', 'low', 'high'])
    
    result = xr.concat([slope, conf_low, conf_high], dim=new_dim)

    # save the result
    sodir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_slope/"
    # create the directory if it does not exist
    if not os.path.exists(sodir):
        os.makedirs(sodir)
    filename = f"plev_50000_{fixed_pattern}_first_JJA_extre_slope.nc"
    result.to_netcdf(sodir + filename)
    return result
# %%
for model in ['MPI_GE','MPI_GE_onepct','CanESM2','CESM1_CAM5','MK36','GFDL_CM3']:
    extrc = read_extrc(model)
    time = '1950-06-01'
    if time != "all" and model != "MPI_GE_onepct":
        time = np.datetime64(time)
        extrc = extrc.sel(time=slice(time, None))
    else:
        extr = extrc
    SMILE_slope(extrc)
# %%
# MPI_GE random
odir = f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_random/extreme_count/"

for ens_size in [20,30,40,50,80]:
    filename = odir + f"plev_50000_decade_JJA_first_{ens_size}_extre_counts.nc"
    ds = xr.open_dataset(filename)
    ds = ds.sel(time = slice('1950-06-01',None))

    # divide the ensemble size 
    ds = ds.pc/ens_size
    sodir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_random/extreme_slope/"
    sfilename = f"plev_50000_decade_JJA_first_{ens_size}_extre_slope.nc"
    slope = SMILE_slope(ds)
    if not os.path.exists(sodir):
        os.makedirs(sodir)
    slope.to_netcdf(sodir + sfilename)
    
# %%
