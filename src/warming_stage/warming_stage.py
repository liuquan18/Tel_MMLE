# %%
import xarray as xr
import numpy as np
import pandas as pd
import warnings
from pandas.tseries.offsets import DateOffset


def read_tsurf_fldmean(dir):
    tsurf = xr.open_dataset(dir,)  # already pre-processed

    try:
        tsurf["time"] = tsurf.indexes["time"].to_datetimeindex()
    except AttributeError:
        pass

    try:
        tsurf = tsurf.tsurf
    except AttributeError:
        tsurf = tsurf.ts

    try:
        tsurf.lon.size == 1 & tsurf.lat.size == 1
    except ValueError:
        print("the fldmean temperature should be calculated first")

    # ens mean
    try:
        fld_ens_mean = tsurf.mean(dim="ens")
    except ValueError:
        fld_ens_mean = tsurf

    # squeeze
    mean = fld_ens_mean.squeeze()
    return mean

# %%

def CO2_period(pc):
    """select the year from pc, 'first10' and 'last10'"""
    years = pc.time
    first10 = slice(years[0], years[9])
    last10 = slice(years[-10], years[-1])
    periods = [first10, last10]
    return periods

def temp_period(fldmean: xr.DataArray):
    """
    warming periods of 0, 2, 4 K, from tsurf fldmean.
    **Arguments**
        *fldmean* the fldmean of tsurf
    **Return**  
        *periods* a list of warming periods
    """
    if isinstance(fldmean, xr.DataArray):
        pass
    else:
        print("only DataArray is accapted, DataSet recevied")

    # anomaly
    period_mean = fldmean.isel(time=slice(0, 10)).mean()  # mean as the basis
    anomaly = fldmean - period_mean
    periods = []

    # 0 degree (1855)
    periods.append(year_to_period(anomaly[5]))
    # 2 degree
    periods.append(
        year_to_period(anomaly.where(anomaly >= 2, drop=True).squeeze()[0])
    )
    # 4 degree
    try:
        periods.append(
            year_to_period(anomaly.where(anomaly >= 4, drop=True).squeeze()[0])
        )
    except IndexError:
        warnings.warn("No fldmean above 4 degree. use the last 10 years instead")
        periods.append(
            slice(
                str(anomaly[-11].time.dt.year.values),
                str(anomaly[-1].time.dt.year.values),
            )
        )
    return periods

def year_to_period(mid_year):
    """return the ten year slice to select"""

    start = mid_year.time.values + DateOffset(years=-4)
    end = mid_year.time.values + DateOffset(years=5)
    return slice(str(start.year), str(end.year))


def split_period(pc,compare, fldmean_tsurf = None):   
    """
    two ways to split the periods, 
    either into 'first10' and 'last10' (compare = 'CO2')
    or into '0', '2', '4' K warming periods (compare = 'temp')
    """
    if compare == "CO2":
        periods = CO2_period(pc)
        period_name = ["first10", "last10"]
    elif (compare == "temp") & (fldmean_tsurf is not None):
        periods = temp_period(fldmean_tsurf)
        period_name = ["0K", "2K", "4K"]
    pcs_period = []
    for i, period in enumerate(periods):
        pc_period = pc.sel(time=period)
        pc_period["compare"] = period_name[i]
        pcs_period.append(pc_period)
    return pcs_period, periods


#
# %%
