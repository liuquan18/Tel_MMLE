#%%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import mpi4py.MPI as MPI
import glob
import os
import sys

# %%

def average_duration(event):
    """
    the average duraion for one single pixel
    """
    
    df = event.to_dataframe()
    df = df['IB index'].reset_index()
    # see https://mp.weixin.qq.com/s?__biz=MzkzNjI1ODQxMQ==&mid=2247484296&idx=1&sn=a9b63336bfbf13724d8b9d13321c75e4&chksm=c2a03d8cf5d7b49afef7ccbf4e0337d019ec1419a241a45ca8e6ca5cfcba6658dc39a9eebd04&token=1605307128&lang=zh_CN#rd
    # Grouper = df.groupby(df.time.dt.year)['IB index'].transform(lambda x: (x<=0).cumsum())
    # G=df[df['IB index']> 0].groupby([df.time.dt.year,Grouper])
    G = df[df['IB index']> 0].groupby((df['IB index']<=0).cumsum())
    durations = G.size()
    dur_mean = durations.mean()
    dur_mean_x = xr.DataArray(dur_mean, dims = ['z'],coords = {'z':event.z})
    return dur_mean_x

def duration_spatial_applyer(event):
    event = event.stack(z=('lat', 'lon'))
    mean_dur = event.groupby('z').apply(average_duration)
    return mean_dur.unstack()

def annual_average_duration(arr):
    annual_dur = arr.resample(time = 'AS').apply(duration_spatial_applyer)
    return annual_dur

def decade_average_duration(arr):
    decade_dur = arr.resample(time = '10AS').apply(duration_spatial_applyer)
    return decade_dur

# %%
def events_count(event, threshold=10):
    """
    count the number of events that lasted more than 10 days for one single pixel
    """
    df = event.to_dataframe()
    df = df['IB index'].reset_index()
    # Grouper = df.groupby(df.time.dt.year)['IB index'].transform(lambda x: (x<=0).cumsum())
    # G=df[df['IB index']> 0].groupby([df.time.dt.year,Grouper])
    G = df[df['IB index']> 0].groupby((df['IB index']<=0).cumsum())
    durations = G.size()
    count = durations[durations> threshold].count()
    count_x = xr.DataArray(count, dims = ['z'],coords = {'z':event.z})
    return count_x

def count_spatial_applyer(arr, threshold=10):
    """
    statistcis (average, count) of the duration of events 
    """
    arr = arr.stack(z=('lat', 'lon'))
    count = arr.groupby('z').apply(events_count, threshold=threshold)
    return count.unstack()

def annual_events_count(arr, threshold=10): # 10 days (6h a day)
    annual_count = arr.resample(time = 'AS').apply(count_spatial_applyer, threshold=threshold)
    return annual_count

def decade_events_count(arr, threshold=40): # 10 days (6h a day)
    decade_count = arr.resample(time = '10AS').apply(count_spatial_applyer, threshold=threshold)
    return decade_count

