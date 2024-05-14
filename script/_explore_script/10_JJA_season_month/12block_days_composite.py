#%%
import xarray as xr
import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import proplot as pplt
import matplotlib as mpl
import src.plots.utils as utils

# %%
# single map for bDays
def bDays_map_single(bDays,ax,levels = np.arange(5,11,1)):
    if levels.min() >= 0:
        cmap = plt.cm.viridis.copy()
        # mask the value of bDays below levels.min()-(levels[1]-levels[0])
        mask_values = levels.min()-2*(levels[1]-levels[0])
        bDays = bDays.where(bDays >= mask_values)
        cmap.set_bad(color = 'none')

    else:
        cmap = plt.cm.RdBu_r.copy()
        interval = levels[1]-levels[0]
        bDays = bDays.where((bDays <=-1*interval) | (bDays >=interval))
        cmap.set_bad(color = 'none')


    p = ax.contourf(
        bDays.lon,
        bDays.lat,
        bDays.values,
        transform=ccrs.PlateCarree(),
        levels=levels,
        cmap=cmap,
        extend = 'both',
    )


    return p, ax

# %%
def bDays_map_composite(bDays,levels_full = np.arange(4,11,1),levels_diff = np.arange(-4,5,1),sig = False):
    try:
        bDays = bDays.sel(mode = 'NAO')
    except:
        pass
    bDays = utils.erase_white_line(bDays)

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
        leftlabels=["pos", "neg"],
        suptitle=f"Change in blocking event occurence in full field of 500 hPa",
        # set the fontsize of labels to 25
    )

    for i, extr_type in enumerate(['pos','neg']):
        Maps = []
        for j, period in enumerate(['first','last','diff']):
            ax = axes[i,j]
            if j<2:
                levels = levels_full
            elif j==2:
                levels = levels_diff
            map, ax = bDays_map_single(bDays.sel(extr_type = extr_type, period = period),
                                        ax = ax,
                                        levels = levels)
            Maps.append(map)
    
            if sig:
                axes[i,2].contourf(
                bDays.sel(extr_type = extr_type, period = 'diff_sig'),
                levels=[-0.5, 0.5, 1.5],
                colors=["none", "none"],
                hatches=["", "xxxx"],
            )

    # Add colorbars for axes[1,1] and axes[1,2]
    cbar1 = fig.colorbar(Maps[1], ax=axes[0,2], orientation='vertical', label='occurence',width = 0.1)
    cbar2 = fig.colorbar(Maps[-1], ax=axes[1,2], orientation='vertical', label='changes in occurence',width = 0.1)
            
    return fig
# %%
bDays_full = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/composite/plev_50000_decade_mpi_first_JJA_JJA_first_last_bDays_composite_mean.nc")
bDays_map_composite(bDays_full['IB index'].sel(mode = 'NAO'),
                    levels_full= np.arange(4,11,1),
                    levels_diff = np.arange(-4,5,1),
                    sig = True)
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/bocu_composite_map_500hPa.png",dpi = 300)
# %%
bDays_ano = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/composite/plev_50000_decade_mpi_first_JJA_JJA_first_last_bDays_ano_composite_mean.nc")
maps = bDays_map_composite(bDays_ano['IB index'].sel(mode = 'NAO'),
                    levels_full= np.arange(0.3,1,0.1),
                    levels_diff = np.arange(-0.4,0.5,0.1),
                    sig = True)
# change the title to anomaly
maps.suptitle("Change in blocking event occurence in anomaly field of 500 hPa")
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/bocu_composite_map_500hPa_remove.png",dpi = 300)

# %%
