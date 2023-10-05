#%%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import proplot as pplt
import cartopy.crs as ccrs
import src.plots.utils as utils
import src.MMLE_TEL.index_stats as index_stats

# %%
def read_wind_months(var_name = 'u', remove_ensmean = False):
    odir = '/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/'
    Jun_fl = odir + f'{var_name}_Jun/'
    Jul_fl = odir + f'{var_name}_Jul/'
    Aug_fl = odir + f'{var_name}_Aug/'

    Jun_f = index_stats.read_var_data(Jun_fl,remove_ensmean=remove_ensmean)
    Jul_f = index_stats.read_var_data(Jul_fl,remove_ensmean=remove_ensmean)
    Aug_f = index_stats.read_var_data(Aug_fl,remove_ensmean=remove_ensmean)

    Jun_f = index_stats.read_var_data(Jun_fl,remove_ensmean=remove_ensmean)
    Jul_f = index_stats.read_var_data(Jul_fl,remove_ensmean=remove_ensmean)
    Aug_f = index_stats.read_var_data(Aug_fl,remove_ensmean=remove_ensmean)

    JJA_f = xr.concat([Jun_f, Jul_f, Aug_f], dim='time')
    JJA_f = JJA_f.sortby('time')

    return JJA_f[var_name]

def wind_map_single(u, v, ax,levels=np.arange(5, 21, 4),cmap = "YlOrRd"):
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
    )
    return contourf, ax


# %%
def duration_map_single(duration, ax, levels = np.arange(5,11,0.5),cmap = 'Fire'):
    
    fig = ax.contourf(
        duration,
        transform=ccrs.PlateCarree(),
        levels=levels,
        cmap = cmap
    )
    return fig, ax

# %%
dec_avedur = xr.open_mfdataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/block_nollb_duration_average/*.nc", combine='nested', concat_dim='ens', parallel=True)
# %%
dec_avedur = dec_avedur['average_duration'].mean(dim = 'ens')
# %%
dec_avedur = utils.erase_white_line(dec_avedur)
# %%
block = [dec_avedur.isel(time = 0)/4, # from 6h to 24h
        dec_avedur.isel(time = -1)/4,
        (dec_avedur.isel(time = -1) - dec_avedur.isel(time = 0))/4]

#%%
u = read_wind_months(var_name = 'u')

v = read_wind_months(var_name = 'v')
u_dec = u.resample(time = '10AS').mean()
v_dec = v.resample(time = '10AS').mean()

u_dec_mean = u_dec.mean(dim = 'ens')
v_dec_mean = v_dec.mean(dim = 'ens')

U = u_dec_mean.sel(plev = 20000).drop('plev')
V = v_dec_mean.sel(plev = 20000).drop('plev')
U= utils.erase_white_line(U)
V = utils.erase_white_line(V)


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
    leftlabels=["wind-200", "blocking"],
    suptitle=f"Change in wind-200 and blocking duration in full field",
    # set the fontsize of labels to 25
)


wind_map, ax = wind_map_single(U.isel(time = 0), V.isel(time = 0), axes[0,0], levels=np.arange(5,31,5))
wind_map_single(U.isel(time = -1), V.isel(time = -1), axes[0,1], levels=np.arange(5,31,5))
wind_map_single(U.isel(time = -1) - U.isel(time = 0), V.isel(time = -1) - V.isel(time = 0), axes[0,2], 
                levels=np.arange(5,31,5))


block_map, ax = duration_map_single(block[0], axes[1,0],levels=np.arange(5,10.5,0.5))
duration_map_single(block[1], axes[1,1],levels=np.arange(5,10.5,0.5))
duration_map_single(block[2], axes[1,2],levels=np.arange(5,10.5,0.5))



fig.colorbar(wind_map, row = 1, loc = 'r', orientation='vertical', label = 'm/s')
fig.colorbar(block_map, row = 2, loc = 'r', orientation='vertical', label = 'days')

fig.save('/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/decadal_wind_blocking.png', dpi=300)
# %%
