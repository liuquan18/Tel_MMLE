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

#%%
def wind_field():
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

    return U, V
