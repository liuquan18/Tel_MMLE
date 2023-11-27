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
def read_data(var_name = 'block_nollb_duration_average'):
    odir = '/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/'
    Jun_fl = odir + f'{var_name}_Jun/'
    Jul_fl = odir + f'{var_name}_Jul/'
    Aug_fl = odir + f'{var_name}_Aug/'

    Jun_f = index_stats.read_var_data(Jun_fl,remove_ensmean=False)
    Jul_f = index_stats.read_var_data(Jul_fl,remove_ensmean=False)
    Aug_f = index_stats.read_var_data(Aug_fl,remove_ensmean=False)

    JJA_f = xr.concat([Jun_f, Jul_f, Aug_f], dim='time')
    JJA_f = JJA_f.sortby('time')
    return JJA_f

#%%
def decadal_change(JJA_f,reduction = 'mean',plev = 20000):
    if reduction == 'mean':
        decade_f = JJA_f.resample(time = '10AS').mean()
        decade_f = decade_f.mean(dim = 'ens') # mean over ensemble
    elif reduction == 'sum':
        decade_f = JJA_f.resample(time = '10AS').sum()
        decade_f = decade_f.sum(dim = 'ens')

    try:
        decade_f = decade_f.sel(plev = plev).drop('plev')
    except:
        pass
    first = decade_f.isel(time = 0)
    last = decade_f.isel(time = -1)
    change = last - first
    return first, last, change

#%%
Duration = read_data(var_name = 'block_nollb_duration_average')
Count = read_data(var_name = 'block_nollb_event_count')

Duration = Duration['duration']
Count = Count['event_count']

Duration = utils.erase_white_line(Duration)
Count = utils.erase_white_line(Count)

#%%
U = read_data(var_name = 'u')
V = read_data(var_name = 'v')

U = U['u']
V = V['v']

U = utils.erase_white_line(U)
V = utils.erase_white_line(V)

#%%
average_durs = decadal_change(Duration,reduction = 'mean')
event_counts = decadal_change(Count,reduction = 'sum')

#%%
average_Us = decadal_change(U,reduction = 'mean',plev = 20000)
average_Vs = decadal_change(V,reduction = 'mean',plev = 20000)


# %%
def block_map_single(duration, ax, levels = np.arange(5,11,0.5),cmap = 'Fire', cut = False):
    if cut:
        cmap = plt.cm.get_cmap(cmap).copy()
        mask_values = levels.min() - (levels[1] - levels[0])
        duration = duration.where(duration >= mask_values)
        cmap.set_bad('none')
    else:
        cmap = cmap
    map = ax.contourf(
        duration.lon,
        duration.lat,
        duration, # from 6h to 1day
        transform=ccrs.PlateCarree(),
        levels=levels,
        cmap = cmap,
    )
    return map, ax

#%%
def wind_map_single(u, v, ax,levels=np.arange(5, 21, 4),cmap = "YlOrRd", color = 'white'):
    # caluclate the wind speed (m/s)
    wind_speed = np.sqrt(u**2 + v**2)

    # plot the wind vector as quiver, and the wind speed as color
    contourf = ax.contourf(
        wind_speed,
        transform=ccrs.PlateCarree(), 
        cmap=cmap, 
        levels = levels,
        extend = 'both',
    )
    ax.quiver(
        u.lon[::4],
        u.lat[::4],
        u.values[::4, ::4],
        v.values[::4, ::4],
        transform=ccrs.PlateCarree(),
        color = color,
    )
    return contourf, ax

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

durationmap,ax = block_map_single(average_durs[0], axes[0,0], 
                               cmap = 'Fire',levels=np.arange(5,10.5,0.5))
block_map_single(average_durs[1], axes[0,1], cmap = 'Fire',levels=np.arange(5,10.5,0.5))
block_map_single(average_durs[2], axes[0,2], cmap = 'Fire',levels=np.arange(5,10.5,0.5))

# first 10 years
wind, ax = wind_plot.wind_map_single(average_Us[0], average_Vs[0], axes[1,0], levels=np.arange(5,31,5),cmap = 'YlOrRd')
countmap, ax = block_map_single(event_counts[0], axes[1,0], cmap = 'viridis',levels=np.arange(40,90.5,10),cut = True)

# last 10 years
wind_plot.wind_map_single(average_Us[1], average_Vs[1], axes[1,1], levels=np.arange(5,31,5),cmap = 'YlOrRd')
block_map_single(event_counts[1], axes[1,1], cmap = 'viridis',levels=np.arange(40,90.5,10),cut = True)

# last 10 - first 10
wind_plot.wind_map_single(average_Us[2],average_Vs[2], axes[1,2],
                            levels=np.arange(5,31,5),cmap = 'YlOrRd')   
block_map_single(event_counts[2], axes[1,2], cmap = 'viridis',levels=np.arange(40,90.5,10),cut=True)


fig.colorbar(durationmap, row = 1, loc = 'r', orientation='vertical', label = 'days',shrink = 0.8)
fig.colorbar(countmap, row = 2, loc = 'r', orientation='vertical', label = 'counts',shrink = 0.8,width = 0.2)
fig.colorbar(wind, row = 2, loc = 'r', orientation='vertical', label = 'm/s',shrink = 0.8,width = 0.2)

# plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/blocking_change.png", dpi = 300)
# %%
