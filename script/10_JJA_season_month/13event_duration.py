#%%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
# %%

def average_duration(arr):
    """
    the average duraion for one single pixel
    """
    
    df = arr.to_dataframe()
    df = df['IB index'].reset_index()
    G=df[df['IB index']> 0].groupby((df['IB index']<=0).cumsum())
    durations = G.size()
    dur_mean = durations.mean()
    dur_mean_x = xr.DataArray(dur_mean, dims = ['z'],coords = {'z':arr.z})
    return dur_mean_x

def yearly_average_duration(arr):
    """
    statistcis (average, count) of the duration of events 
    """
    arr = arr.stack(z=('lat', 'lon'))
    mean_dur = arr.groupby('z').apply(average_duration)
    return mean_dur.unstack()

def decade_average_duration(arr):
    decade_dur = arr.resample(time = '10AS').apply(yearly_average_duration)
    return decade_dur.mean(dim = 'ens')

# %%
def events_count(arr, threshold=10):
    """
    count the number of events that lasted more than 10 days for one single pixel
    """
    df = arr.to_dataframe()
    df = df['IB index'].reset_index()
    G=df[df['IB index']> 0].groupby((df['IB index']<=0).cumsum())
    durations = G.size()
    count = durations[durations> threshold].count()
    count_x = xr.DataArray(count, dims = ['z'],coords = {'z':arr.z})
    return count_x

def yearly_events_count(arr, threshold=10):
    """
    statistcis (average, count) of the duration of events 
    """
    arr = arr.stack(z=('lat', 'lon'))
    count = arr.groupby('z').apply(events_count, threshold=threshold)
    return count.unstack()

def decade_events_count(arr, threshold=10):
    decade_count = arr.resample(time = '10AS').apply(yearly_events_count, threshold=threshold)
    return decade_count.sum(dim = 'ens')

#%%
events = xr.open_mfdataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_event_ano/*.nc", combine='nested', concat_dim='ens')
average_dur = decade_average_duration(events)
average_dur.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/block_event/average_dur_ano.nc")
# %%
events = xr.open_mfdataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_event_ano/*.nc", combine='nested', concat_dim='ens')
event_count = decade_events_count(events)
event_count.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/block_event/event_count_ano.nc")
# %%
#%%
events = xr.open_mfdataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_event/*.nc", combine='nested', concat_dim='ens')
average_dur = decade_average_duration(events)
average_dur.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/block_event/average_dur.nc")
# %%
events = xr.open_mfdataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_event/*.nc", combine='nested', concat_dim='ens')
event_count = decade_events_count(events)
event_count.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/block_event/event_count.nc")
# %%
