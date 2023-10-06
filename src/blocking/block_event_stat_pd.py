#%%
import xarray as xr
import pandas as pd


# %%
def annual_blocking_stat(ix, threshod=40):
    """
    seperate the year, lon, lat from the index of the xarray
    """
    # find the blocks
    df = ix.to_dataframe().reset_index()
    Grouper = df.groupby([df.time.dt.year, 'lat','lon'])['IB index'].transform(lambda x: (x<=0).cumsum())
    G = df[df['IB index']>0].groupby([df.time.dt.year,'lat','lon'])['IB index']
    blocks = G.size()
    blocks = blocks.reset_index()

    # duration
    duration = blocks.groupby(['time','lat','lon'])['IB index'].apply(average_duration)
    Duration = duration.to_xarray()
    Duration.name = 'duration'

    # count of events that lasted more than 10 days
    counts = blocks.groupby(['time','lat','lon'])['IB index'].apply(event_count,threshod=threshod)
    Count = counts.to_xarray()
    Count.name = 'event_count'

    return Duration, Count

def event_count(duration, threshod = 40):
    count = duration[duration>threshod].count()
    return count

def average_duration(blocks):
    return blocks.mean()