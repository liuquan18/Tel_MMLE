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

# %%
zg_first = xr.open_dataarray(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_first_mean.nc"
)
zg_last = xr.open_dataarray(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_last_mean.nc"
)
# %%
zg_diff = zg_last - zg_first

# %%
zg_first = erase_white_line(zg_first)
zg_last = erase_white_line(zg_last)
zg_diff = erase_white_line(zg_diff)
# %%
zg_first = zg_first / 1000  # from m to km
zg_last = zg_last / 1000
zg_diff = zg_diff / 1000
# %%
u_first = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/u_1850-1859.nc"
)
u_last = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/u_2090-2099.nc"
)
# %%
u_first = u_first.var131.squeeze()
u_last = u_last.var131.squeeze()
u_diff = u_last - u_first
# %%
u_first = erase_white_line(u_first)
u_last = erase_white_line(u_last)
u_diff = erase_white_line(u_diff)

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
eof_result_rm = xr.open_dataset(
        "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result.nc"
    )

eof_first_rm = eof_result_rm.eof.isel(decade=0).sel(mode = 'NAO')
eof_last_rm = eof_result_rm.eof.isel(decade=-1).sel(mode = 'NAO')

eof_first_rm = erase_white_line(eof_first_rm)
eof_last_rm = erase_white_line(eof_last_rm)
eof_diff_rm = eof_last_rm - eof_first_rm
# %%
fig, axes = plt.subplots(
    nrows=4,
    ncols=3,
    figsize=(12, 12),
    subplot_kw={"projection": ccrs.Orthographic(-20, 60)},
)
zg_levels_seq = np.arange(5, 6, 0.1)
u_levels_seq = np.arange(0, 60, 5)
eof_levels_seq = np.arange(38, 41.5, 0.5)

zg_levels_div = np.arange(-0.2, 0.21, 0.02)
u_levels_div = np.arange(-2, 2.1, 0.5)
eof_levels_div = np.arange(-30, 31, 5)
# Plot zg_first
cf = axes[0, 0].contourf(
    zg_first.lon, zg_first.lat, zg_first, 
    cmap="viridis", 
    levels=zg_levels_seq,
    extend='both',
    transform=ccrs.PlateCarree()
)
axes[0, 0].set_title("first", fontsize=15)
fig.colorbar(cf, ax=axes[0, 0], orientation='vertical', pad=0.05, label="Z500 [km]")

# Plot zg_last
cf = axes[0, 1].contourf(
    zg_last.lon, zg_last.lat, zg_last, 
    levels=zg_levels_seq,
    cmap="viridis", 
    extend='both',
    transform=ccrs.PlateCarree()
)
axes[0, 1].set_title("last", fontsize=15)
fig.colorbar(cf, ax=axes[0, 1], orientation='vertical', pad=0.05, label="Z500 [km]")

# Plot zg_diff
cf = axes[0, 2].contourf(
    zg_diff.lon, zg_diff.lat, zg_diff, 
    levels=zg_levels_div,
    cmap="coolwarm", 
    extend='both',
    transform=ccrs.PlateCarree()
)
axes[0, 2].set_title("last - first", fontsize=15)
fig.colorbar(cf, ax=axes[0, 2], orientation='vertical', pad=0.05, label="Z500 [km]")

# u as contour
axes[0, 0].contour(
    u_first.lon,
    u_first.lat,
    u_first,
    colors="black",
    levels=[level for level in u_levels_seq if level != 0],
    transform=ccrs.PlateCarree(),
)
axes[0, 1].contour(
    u_last.lon,
    u_last.lat,
    u_last,
    colors="black",
    levels=[level for level in u_levels_seq if level != 0],
    transform=ccrs.PlateCarree(),
)
axes[0, 2].contour(
    u_diff.lon,
    u_diff.lat,
    u_diff,
    colors="black",
    levels=[level for level in u_levels_div if level != 0],
    transform=ccrs.PlateCarree(),
)

# eof1 in the second row
cf = axes[1, 0].contourf(
    eof_first.lon,
    eof_first.lat,
    eof_first.isel(mode=0),
    cmap="Reds",
    levels=np.arange(35.5, 39, 0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[1, 0], orientation='vertical', pad=0.05, label="EOF1")


cf  = axes[1, 1].contourf(
    eof_last.lon,
    eof_last.lat,
    eof_last.isel(mode=0),
    cmap="Reds",
    levels=np.arange(35.5, 39, 0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[1, 1], orientation='vertical', pad=0.05, label="EOF1")

cf = axes[1, 2].contourf(
    eof_diff.lon,
    eof_diff.lat,
    eof_diff.isel(mode=0),
    cmap="coolwarm",
    levels=np.arange(-2, 2.1, 0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[1, 2], orientation='vertical', pad=0.05, label="EOF1")

# eof2 in the third row
cf = axes[2, 0].contourf(
    eof_first.lon,
    eof_first.lat,
    eof_first.isel(mode=1),
    cmap="coolwarm",
    levels=np.arange(-2,2.1,0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[2, 0], orientation='vertical', pad=0.05, label="EOF2")

cf = axes[2, 1].contourf(
    eof_last.lon,
    eof_last.lat,
    eof_last.isel(mode=1),
    cmap="coolwarm",
    levels=np.arange(-2,2.1,0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)

fig.colorbar(cf, ax=axes[2, 1], orientation='vertical', pad=0.05, label="EOF2")


cf = axes[2, 2].contourf(
    eof_diff.lon,
    eof_diff.lat,
    eof_diff.isel(mode=1),
    cmap="coolwarm",
    levels=np.arange(-2,2.1,0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[2, 2], orientation='vertical', pad=0.05, label="EOF2")

# forth row for eof_rm
cf = axes[3, 0].contourf(
    eof_first_rm.lon,
    eof_first_rm.lat,
    eof_first_rm,
    cmap="coolwarm",
    levels= np.arange(-2,2.1,0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)

fig.colorbar(cf, ax=axes[3, 0], orientation='vertical', pad=0.05, label="EOF_NAO")

cf = axes[3, 1].contourf(
    eof_last_rm.lon,
    eof_last_rm.lat,
    eof_last_rm,
    cmap="coolwarm",
    levels= np.arange(-2,2.1,0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)

fig.colorbar(cf, ax=axes[3, 1], orientation='vertical', pad=0.05, label="EOF_NAO")


cf = axes[3, 2].contourf(
    eof_diff_rm.lon,
    eof_diff_rm.lat,
    eof_diff_rm,
    cmap="coolwarm",
    levels= np.arange(-2,2.1,0.5),
    extend='both',
    transform=ccrs.PlateCarree(),
)
fig.colorbar(cf, ax=axes[3, 2], orientation='vertical', pad=0.05, label="EOF_NAO")




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
axes[2, 0].text(-0.1, 1.05, "g", transform=axes[2, 0].transAxes, fontsize=15, fontweight='bold')
axes[2, 1].text(-0.1, 1.05, "h", transform=axes[2, 1].transAxes, fontsize=15, fontweight='bold')
axes[2, 2].text(-0.1, 1.05, "i", transform=axes[2, 2].transAxes, fontsize=15, fontweight='bold')
axes[3, 0].text(-0.1, 1.05, "j", transform=axes[3, 0].transAxes, fontsize=15, fontweight='bold')
axes[3, 1].text(-0.1, 1.05, "k", transform=axes[3, 1].transAxes, fontsize=15, fontweight='bold')
axes[3, 2].text(-0.1, 1.05, "l", transform=axes[3, 2].transAxes, fontsize=15, fontweight='bold')


plt.tight_layout()
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/pattern_disscusion.pdf")
# %%
