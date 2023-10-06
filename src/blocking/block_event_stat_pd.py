#%%
import xarray as xr
import pandas as pd


# %%
def annual_blocking_stat(ix, month,threshod=40):
    """
    seperate the year, lon, lat from the index of the xarray
    *month* for chaning the time
    """
    # find the blocks
    df = ix.to_dataframe().reset_index()
    Grouper = df.groupby([df.time.dt.year, 'lat','lon'])['IB index'].transform(lambda x: (x<=0).cumsum())
    G = df[df['IB index']>0].groupby([df.time.dt.year,'lat','lon'])['IB index']
    blocks = G.size()
    blocks = blocks.reset_index()
    dates = pd.to_datetime([year + '-' + month + '-01' for year in blocks.time.values.astype(str)])
    blocks['time'] = dates

    # duration
    duration = blocks.groupby([blocks.time.dt.year,'lat','lon'])['IB index'].apply(average_duration)
    Duration = duration.to_xarray()
    Duration.name = 'duration'
    Duration['time'] = pd.to_datetime([year + '-' + month + '-01' for year in Duration.time.values.astype(str)])

    # count of events that lasted more than 10 days
    counts = blocks.groupby([blocks.time.dt.year,'lat','lon'])['IB index'].apply(event_count,threshod=threshod)
    Count = counts.to_xarray()
    Count.name = 'event_count'
    Count['time'] = pd.to_datetime([year + '-' + month + '-01' for year in Count.time.values.astype(str)])

    return Duration, Count

def event_count(duration, threshod = 40):
    count = duration[duration>threshod].count()
    return count

def average_duration(blocks):
    return blocks.mean()