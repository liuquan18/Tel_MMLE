# %%
import xarray as xr
import numpy as np
import pandas as pd
import statsmodels.api as sm
import os


# %% calculate the slope of the extreme counts
def cal_slope(extrc, tsurf=None):
    if tsurf is not None:
        x = tsurf
        x = sm.add_constant(x)
    else:
        x = np.arange(len(extrc))
        x = sm.add_constant(x)
    y = extrc
    model = sm.OLS(y, x)
    results = model.fit()
    slope = results.params[1]
    return slope


def cal_conf_low(extrc, tsurf=None):
    # calculate the confidence interval of the slope
    if tsurf is not None:
        x = tsurf
        x = sm.add_constant(x)
    else:
        x = np.arange(len(extrc))
        x = sm.add_constant(x)
    y = extrc
    model = sm.OLS(y, x)
    results = model.fit()
    confint = results.conf_int(alpha=0.05, cols=None)[1][0]  # the lower bound
    return confint


def cal_conf_high(extrc, tsurf=None):
    if tsurf is not None:
        x = tsurf
        x = sm.add_constant(x)
    else:
        x = np.arange(len(extrc))
        x = sm.add_constant(x)
    y = extrc
    model = sm.OLS(y, x)
    results = model.fit()
    confint = results.conf_int(alpha=0.05, cols=None)[1][1]  # the upper bound
    return confint


# %%
def read_extrc(model, fixed_pattern="decade_mpi"):
    """read extreme counts"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_count/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_extre_counts.nc"
    ds = xr.open_dataset(odir + filename).pc

    # divide the ensemble size of each model
    ens_sizes = {
        "MPI_GE": 100,
        "MPI_GE_onepct": 100,
        "CanESM2": 50,
        "CESM1_CAM5": 40,
        "MK36": 30,
        "GFDL_CM3": 20,
    }
    ds = ds / ens_sizes[model]
    return ds


# %%
def read_tsurf(model):
    fname = (
        "/work/mh0033/m300883/Tel_MMLE/data/"
        + model
        + "/extreme_count/ens_fld_year_mean.nc"
    )
    ds = xr.open_dataset(fname)
    try:
        tsurf = ds.tsurf
    except AttributeError:
        tsurf = ds.ts
    return tsurf


# %%
def SMILE_slope(extrc, tsurf=None, fixed_pattern="decade_mpi"):
    try:
        extrc = extrc.sel(confidence="true").drop("confidence")
    except KeyError:
        pass
    if tsurf is not None:
        input_dim = [["time"], ["time"]]
    else:
        input_dim = [["time"], []]
    slope = xr.apply_ufunc(
        cal_slope,
        extrc,
        tsurf,
        input_core_dims=input_dim,
        vectorize=True,
        dask="parallelized",
    )
    conf_low = xr.apply_ufunc(
        cal_conf_low,
        extrc,
        tsurf,
        input_core_dims=input_dim,
        vectorize=True,
        dask="parallelized",
    )
    conf_high = xr.apply_ufunc(
        cal_conf_high,
        extrc,
        tsurf,
        input_core_dims=input_dim,
        vectorize=True,
        dask="parallelized",
    )
    new_dim = xr.IndexVariable("slopes", ["true", "low", "high"])

    result = xr.concat([slope, conf_low, conf_high], dim=new_dim)
    return result


# %%
def rate_increase(fixed_pattern = "decade_mpi",x = 'tsurf',time = "1950-06-01"):
    for model in ["MPI_GE", "MPI_GE_onepct", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]:
        if x == 'tsurf':
            tsurf = read_tsurf(model)
        else:
            tsurf = None
        extrc = read_extrc(model)
        if time != "all" and model != "MPI_GE_onepct":
            time = np.datetime64(time)
            extrc = extrc.sel(time=slice(time, None))
            try:
                tsurf = tsurf.sel(time=slice(time, None))
            except:
                pass
        else:
            extrc = extrc
        try:
            tsurf['time'] = extrc['time']
        except: # None tsurf
            pass
        slopes = SMILE_slope(extrc,tsurf=tsurf)

    # save the result
        sodir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_slope/"
    # create the directory if it does not exist
        if not os.path.exists(sodir):
            os.makedirs(sodir)
        filename = f"plev_50000_{fixed_pattern}_first_JJA_extre_slope_{x}.nc"
        slopes.to_netcdf(sodir + filename)

#%%
rate_increase(x = 'tsurf',time = "1950-06-01")
#%%
rate_increase(x = 'time',time = "1950-06-01")

# %%
# MPI_GE random
def resample_rate_increase(x = 'tsurf',time = "1950-06-01"):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_random/extreme_count/"

    for ens_size in [20, 30, 40, 50, 80]:
        filename = odir + f"plev_50000_decade_JJA_first_{ens_size}_extre_counts.nc"
        ds = xr.open_dataset(filename)
        start_time = np.datetime64(time)
        ds = ds.sel(time=slice(start_time, None))
    if x == 'tsurf':
        tsurf = read_tsurf('MPI_GE')
        tsurf = tsurf.sel(time=slice(start_time, None))
        tsurf['time'] = ds['time']
    else:
        tsurf = None

    # divide the ensemble size
        ds = ds.pc / ens_size
        sodir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_random/extreme_slope/"
        sfilename = f"plev_50000_decade_JJA_first_{ens_size}_extre_slope_{x}.nc"
        slope = SMILE_slope(ds,tsurf)
        if not os.path.exists(sodir):
            os.makedirs(sodir)
        slope.to_netcdf(sodir + sfilename)
#%%
resample_rate_increase(x = 'tsurf',time = "1950-06-01")
#%%
resample_rate_increase(x = 'time',time = "1950-06-01")
# %%
