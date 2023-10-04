#%%
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import proplot as pplt
import cartopy.crs as ccrs
import src.plots.utils as utils

# %%
dec_avedur = xr.open_mfdataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30/block_duration_average/*.nc", combine='nested', concat_dim='ens', parallel=True)
# %%
dec_avedur = dec_avedur['average_duration'].mean(dim = 'ens')
# %%
dec_avedur = utils.erase_white_line(dec_avedur)
# %%
data = [dec_avedur.isel(time = 0),
        dec_avedur.isel(time = -1),
        dec_avedur.isel(time = -1) - dec_avedur.isel(time = 0)]
# %%
def duration_map_single(duration, ax, levels = np.arange(5,26,5)):
    pass
# %%
fig = pplt.figure(figsize=(50 / 25.4, 50 / 25.4),sharex=False,sharey=False)
fig.format(
    abcloc="ul",
    abc="a",
)
ax = fig.subplots(
    ncols=1,
    nrows=1,
    proj="ortho",
    proj_kw=({"lon_0": -20, "lat_0": 60})
)

ax.format(
    latlines=20,
    lonlines=30,
    color="grey7",
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
)

duration = dec_avedur.isel(time = 0)
levels = np.arange(5,26,2)
ax.contourf(
    duration,
    transform=ccrs.PlateCarree(),
    levels=levels,
    extend = 'both'
)
# %%
