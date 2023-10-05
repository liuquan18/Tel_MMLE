#%%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import proplot as pplt
import cartopy.crs as ccrs
import src.plots.utils as utils
import src.MMLE_TEL.index_stats as index_stats
import src.plots.wind_plot as wind_plot
# %%
import importlib
importlib.reload(wind_plot)
# %%
def block_map_single(duration, ax, levels = np.arange(5,11,0.5),cmap = 'Fire'):
    
    map = ax.contourf(
        duration,
        transform=ccrs.PlateCarree(),
        levels=levels,
        cmap = cmap,
    )
    return map, ax
# %%
event_count = xr.open_mfdataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/block_nollb_event_count/*.nc",combine='nested', concat_dim='ens', parallel=True)
# %%
average_dur = xr.open_mfdataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/block_nollb_duration_average/*.nc", combine='nested', concat_dim='ens', parallel=True)

#%%
event_count = event_count.sum(dim = 'ens')
average_dur = average_dur.mean(dim = 'ens')
#%%
event_count = utils.erase_white_line(event_count['block_event_count'])
average_dur = utils.erase_white_line(average_dur['average_duration'])

#%%
U,V = wind_plot.wind_field()

# %%
event_counts = [
    event_count.isel(time = 0),
    event_count.isel(time = -1),
    event_count.isel(time = -1) - event_count.isel(time = 0)
]

average_durs = [
    average_dur.isel(time = 0)/4, # day
    average_dur.isel(time = -1)/4, # day
    (average_dur.isel(time = -1) - average_dur.isel(time = 0))/4 # day
]


# %%
fig = pplt.figure(figsize=(150 / 25.4, 100 / 25.4),sharex=False,sharey=False)
fig.format(
    abcloc="ul",
    abc="a",
)
axes = fig.subplots(
    ncols=3,
    nrows=2,
    proj="ortho",
    proj_kw=({"lon_0": -20, "lat_0": 60})
)

axes.format(
    latlines=20,
    lonlines=30,
    color="grey7",
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
    toplabels=["first10", "last10", "last10 - first10"],
    leftlabels=["average duration day", "count of >10 days-event"],
    suptitle=f"Change in blocking in full field",
    # set the fontsize of labels to 25
)

duration,ax = block_map_single(average_durs[0], axes[0,0], 
                               cmap = 'Fire',levels=np.arange(5,10.5,0.5))
block_map_single(average_durs[1], axes[0,1], cmap = 'Fire',levels=np.arange(5,10.5,0.5))
block_map_single(average_durs[2], axes[0,2], cmap = 'Fire',levels=np.arange(5,10.5,0.5))

wind, ax = wind_plot.wind_map_single(U.isel(time = 0), V.isel(time = 0), axes[1,0], levels=np.arange(5,31,5),cmap = 'YlOrRd')
count,ax = block_map_single(event_counts[0], axes[1,0], 
                            cmap = 'viridis',levels=np.arange(20,60.5,5))

wind_plot.wind_map_single(U.isel(time = -1), V.isel(time = -1), axes[1,1], levels=np.arange(5,31,5),cmap = 'YlOrRd')
block_map_single(event_counts[1], axes[1,1], cmap = 'viridis',levels=np.arange(20,60.5,5))

wind_plot.wind_map_single(U.isel(time = -1) - U.isel(time = 0), V.isel(time = -1) - V.isel(time = 0), axes[1,2],
                            levels=np.arange(5,31,5),cmap = 'YlOrRd')   
block_map_single(event_counts[2], axes[1,2], cmap = 'viridis',levels=np.arange(20,60.5,5))


fig.colorbar(duration, row = 1, loc = 'r', orientation='vertical', label = 'days',shrink = 0.8)
fig.colorbar(count, row = 2, loc = 'r', orientation='vertical', label = 'counts',shrink = 0.8,width = 0.2)
fig.colorbar(wind, row = 2, loc = 'r', orientation='vertical', label = 'm/s',shrink = 0.8,width = 0.2)

# plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/blocking_change.png", dpi = 300)
# %%
