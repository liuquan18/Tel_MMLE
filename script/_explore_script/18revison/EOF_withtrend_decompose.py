# %%
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
eof_result_dec = xr.open_dataset(
        "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result_withtrend.nc"
    )

eof_first = eof_result_dec.eof.isel(decade=0)

eof_last = eof_result_dec.eof.isel(decade=-1)

eof_first = -1 * erase_white_line(eof_first)
eof_last = -1 * erase_white_line(eof_last)

eof_diff = eof_last - eof_first
#%%
eof_result_all = xr.open_dataset(
        "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/plev_50000_all_first_JJA_eof_result_withtrend.nc"
    )

eof_all_one = eof_result_all.eof.isel(mode=0, decade = 0)
eof_all_two = eof_result_all.eof.isel(mode=1, decade = 0)

eof_all_one = -1 * erase_white_line(eof_all_one)
eof_all_two = -1 * erase_white_line(eof_all_two)
# %%
fig, axes = plt.subplots(
    nrows=2,
    ncols=3,
    figsize=(12, 8),
    subplot_kw={"projection": ccrs.Orthographic(-20, 60)},
)
zg_levels_seq = np.arange(5, 6, 0.1)
u_levels_seq = np.arange(0, 60, 5)
eof_levels_seq = np.arange(38, 41.5, 0.5)

zg_levels_div = np.arange(-0.2, 0.21, 0.02)
u_levels_div = np.arange(-2, 2.1, 0.5)
eof_levels_div = np.arange(-30, 31, 5)

# eof1 in the first row
cf = axes[0, 0].contourf(
    eof_first.lon,
    eof_first.lat,
    eof_first.isel(mode=0),
    cmap="Reds",
    levels=np.arange(35.5, 39, 0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[0, 0], orientation='vertical', pad=0.05, label="EOF1")


cf  = axes[0, 1].contourf(
    eof_last.lon,
    eof_last.lat,
    eof_last.isel(mode=0),
    cmap="Reds",
    levels=np.arange(35.5, 39, 0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[0, 1], orientation='vertical', pad=0.05, label="EOF1")

cf = axes[0, 2].contourf(
    eof_diff.lon,
    eof_diff.lat,
    eof_diff.isel(mode=0),
    cmap="coolwarm",
    levels=np.arange(-2, 2.1, 0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[0, 2], orientation='vertical', pad=0.05, label="EOF1")

# eof2 in the second row
cf = axes[1, 0].contourf(
    eof_first.lon,
    eof_first.lat,
    eof_first.isel(mode=1),
    cmap="coolwarm",
    levels=np.arange(-2,2.1,0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[1, 0], orientation='vertical', pad=0.05, label="EOF2")

cf = axes[1, 1].contourf(
    eof_last.lon,
    eof_last.lat,
    eof_last.isel(mode=1),
    cmap="coolwarm",
    levels=np.arange(-2,2.1,0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)

fig.colorbar(cf, ax=axes[1, 1], orientation='vertical', pad=0.05, label="EOF2")


cf = axes[1, 2].contourf(
    eof_diff.lon,
    eof_diff.lat,
    eof_diff.isel(mode=1),
    cmap="coolwarm",
    levels=np.arange(-2,2.1,0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[1, 2], orientation='vertical', pad=0.05, label="EOF2")

for ax in axes.flat:
    ax.coastlines()
    ax.set_global()

# add a, b, c
axes[0, 0].text(-0.1, 1.05, "a", transform=axes[0, 0].transAxes, fontsize=15, fontweight='bold')
axes[0, 1].text(-0.1, 1.05, "b", transform=axes[0, 1].transAxes, fontsize=15, fontweight='bold')
axes[0, 2].text(-0.1, 1.05, "c", transform=axes[0, 2].transAxes, fontsize=15, fontweight='bold')
axes[1, 0].text(-0.1, 1.05, "d", transform=axes[1, 0].transAxes, fontsize=15, fontweight='bold')
axes[1, 1].text(-0.1, 1.05, "e", transform=axes[1, 1].transAxes, fontsize=15, fontweight='bold')
axes[1, 2].text(-0.1, 1.05, "f", transform=axes[1, 2].transAxes, fontsize=15, fontweight='bold')


plt.tight_layout()
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/EOF_withtrend.pdf")
# %%
