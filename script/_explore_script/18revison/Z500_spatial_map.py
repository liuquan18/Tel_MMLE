#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import src.MMLE_TEL.index_stats as index_stats
import glob
from src.plots.utils import erase_white_line
import proplot as pplt
# %%
# odir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg/"
# # %%
# zg = xr.open_mfdataset(odir + "*.nc", combine='nested', concat_dim='ens', parallel=True)
# # %%
# zg_first = zg.sel(time = slice('1850','1859'), plev = 50000).var156
# # %%
# zg_last = zg.sel(time = slice('2090', '2099'), plev = 50000).var156
# # %%
# zg_first_mean = zg_first.mean(dim=('time','ens'))
# zg_last_mean = zg_last.mean(dim=('time','ens'))

# # %%
# zg_first_mean.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_first_mean.nc")
# zg_last_mean.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_last_mean.nc")
#%%
zg_first = xr.open_dataarray("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_first_mean.nc")
zg_last = xr.open_dataarray("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_last_mean.nc")
#%%
zg_diff = zg_last - zg_first

#%%
zg_first = erase_white_line(zg_first)
zg_last = erase_white_line(zg_last)
zg_diff = erase_white_line(zg_diff)
#%%
zg_first = zg_first/1000 #from m to km
zg_last = zg_last/1000
zg_diff = zg_diff/1000

# %%
fig3, axes = pplt.subplots(
    width=180 / 25.4,
    proj="ortho",
    proj_kw=({"lon_0": -20, "lat_0": 60}),
    nrows=1,
    ncols=3,
)
axes.format(
    abc = True,
    latlines=20,
    lonlines=30,
    color="grey7",
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
    toplabels=["first", "last", "last - first"],
    toplabels_kw={"fontsize": 7, },
)

zg_first.plot.contourf(
    ax=axes[0],
    transform=ccrs.PlateCarree(),
    cmap="viridis",
    levels=np.arange(5,6,0.1),
    extend="both",
    cbar=True,
    cbar_kw={"label": "Z500 [km]"},
)
zg_last.plot.contourf(
    ax=axes[1],
    transform=ccrs.PlateCarree(),
    cmap="viridis",
    levels=np.arange(5,6,0.1),
    extend="both",
    cbar=True,
    cbar_kw={"label": "Z500 [km]"},
)

zg_diff.plot.contourf(
    ax=axes[2],
    transform=ccrs.PlateCarree(),
    cmap="coolwarm",
    levels=np.arange(-0.2,0.21,0.02),
    extend="both",
    cbar=True,
    cbar_kw={"label": "Z500 [km]"},
)



# %%
