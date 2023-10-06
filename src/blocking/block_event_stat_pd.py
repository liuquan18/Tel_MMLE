#%%
import xarray as xr
import pandas as pd


# %%
def annual_blocking_stat(ix, threshod=40):
    """
    seperate the year, lon, lat from the index of the xarray
    """
    df = ix.to_dataframe().reset_index()
    Grouper = df.groupby([df.time.dt.year, 'lat','lon'])['IB index'].transform(lambda x: (x<=0).cumsum())
    G = df[df['IB index']>0].groupby([df.time.dt.year,'lat','lon'])['IB index']

    blocks = G.size()
    # duration
    duration = G.size()
    Duration = duration.to_xarray()

    # count of events that lasted more than 10 days
    duration = duration.reset_index()
    duration = duration.rename(columns = {'IB index':'duration'})
    counts = duration.groupby(['time','lat','lon'])['duration'].apply(count,threshod=threshod)
    Count = counts.to_xarray()

    return Duration, Count

def event_count(duration, threshod = 40):
    count = duration[duration>threshod].count()
    return count

def average_duration(blocks):
    return blocks.mean()