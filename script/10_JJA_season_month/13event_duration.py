#%%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
# %%
event = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_event_ano/onepct_1850-1999_ens_0069.gph500.nc")
event = event['IB index']
# %%
first = event.sel(time = event.time.dt.year == 1950)
# %%
df = ex.to_dataframe().reset_index()
# %%
df = df[['time','IB index']]
# %%
G=df[df['IB index']> 0].groupby((df['IB index']<=0).cumsum())
# %%
Events = G.agg(
   start_time = pd.NamedAgg(column = 'time',aggfunc='min'),
   duration = pd.NamedAgg(column = 'time',aggfunc = 'size'),
   sum = pd.NamedAgg(column = 'IB index',aggfunc = 'sum'),)
# %%

def average_duration(arr):
    """
    the average duraion for one single pixel
    """
    
    df = arr.to_dataframe().reset_index()
    G=df[df['IB index']> 0].groupby((df['IB index']<=0).cumsum())
    durations = G.size()
    return durations.mean()
# %%
def events_count(arr, threshold=10):
    """
    count the number of events that lasted more than 10 days for one single pixel
    """
    df = arr.to_dataframe().reset_index()
    G=df[df['IB index']> 0].groupby((df['IB index']<=0).cumsum())
    durations = G.size()
    return durations[durations> threshold].count()
# %%
# apply the function to the whole dataset

first_average_dur = xr.apply_ufunc(average_duration, 
                                   first, 
                                   input_core_dims=[['time']],
                                   output_core_dims=
                                   dask='parallelized', 
                                   output_dtypes=[float])
# %%
