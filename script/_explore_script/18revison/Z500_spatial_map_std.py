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
# zg = xr.open_mfdataset(odir + "*.nc", combine='nested', concat_dim='ens', parallel=True)

# zg_first = zg.sel(time = slice('1850','1859'), plev = 50000).var156
# zg_last = zg.sel(time = slice('2090', '2099'), plev = 50000).var156

# zg_first_std = zg_first.std(dim=('time','ens'))
# zg_last_std = zg_last.std(dim=('time','ens'))

# zg_first_std.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_first_std.nc")
# zg_last_std.to_netcdf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_last_std.nc")
# %%
zg_first_std = xr.open_dataarray(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_first_std.nc"
)
zg_last_std = xr.open_dataarray(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_last_std.nc"
)
zg_diff_std = zg_last_std - zg_first_std

# %%
zg_first_std = erase_white_line(zg_first_std)
zg_last_std = erase_white_line(zg_last_std)
zg_diff_std = erase_white_line(zg_diff_std)
# %%
# %%
zg_first_mean = xr.open_dataarray(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_first_mean.nc"
)
zg_last_mean = xr.open_dataarray(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_first_last/zg_last_mean.nc"
)
zg_diff_mean = zg_last_mean - zg_first_mean

# %%
zg_first_mean = erase_white_line(zg_first_mean)
zg_last_mean = erase_white_line(zg_last_mean)
zg_diff_mean = erase_white_line(zg_diff_mean)
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
u_first_nozon = u_first - u_first.mean(dim = 'lon')
u_last_nozon = u_last - u_last.mean(dim = 'lon')
u_diff_nozon = u_diff - u_diff.mean(dim = 'lon')


#%%
# the projected spatial pattern
eof_first_rm = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/first_pattern_projected.nc"
).__xarray_dataarray_variable__
eof_last_rm = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/last_pattern_projected.nc"
).__xarray_dataarray_variable__

eof_first_rm = erase_white_line(eof_first_rm)
eof_last_rm = erase_white_line(eof_last_rm)
eof_diff_rm = eof_last_rm - eof_first_rm

# %%
fig, axes = plt.subplots(
    nrows=3,
    ncols=3,
    figsize=(12, 12),
    subplot_kw={"projection": ccrs.Orthographic(-20, 60)},
)

zg_levels_mean = np.arange(5000, 6000, 100)
zg_levels_mean_div = np.arange(-180, 181, 30)

zg_levels_seq = np.arange(50, 200, 10)
u_levels_seq = np.arange(0, 60, 5)
u_levels_seq_nonzon = np.arange(-4.5, 5, 1.5)

zg_levels_div = np.arange(-30, 31, 5)
u_levels_div = np.arange(-1.5, 1.6, 0.5)
eof_levels_div = np.arange(-30, 31, 5)
# First row for zg_mean
# Plot zg_first_mean
cf = axes[0, 0].contourf(
    zg_first_mean.lon, zg_first_mean.lat, zg_first_mean, 
    cmap="viridis", 
    levels=zg_levels_mean,
    extend='both',
    transform=ccrs.PlateCarree()
)
axes[0, 0].set_title("first10 Z500 mean", fontsize=15)
fig.colorbar(cf, ax=axes[0, 0], orientation='vertical', pad=0.05, label="Z500 [m]", shrink=0.8)

# Plot zg_last_mean
cf = axes[0, 1].contourf(
    zg_last_mean.lon, zg_last_mean.lat, zg_last_mean, 
    levels=zg_levels_mean,
    cmap="viridis", 
    extend='both',
    transform=ccrs.PlateCarree()
)
axes[0, 1].set_title("last10 Z500 mean", fontsize=15)
fig.colorbar(cf, ax=axes[0, 1], orientation='vertical', pad=0.05, label="Z500 [m]", shrink=0.8)

# Plot zg_diff_mean
cf = axes[0, 2].contourf(
    zg_diff_mean.lon, zg_diff_mean.lat, zg_diff_mean, 
    levels=zg_levels_mean_div,
    cmap="coolwarm", 
    extend='both',
    transform=ccrs.PlateCarree()
)
axes[0, 2].set_title("diff Z500 mean", fontsize=15)
fig.colorbar(cf, ax=axes[0, 2], orientation='vertical', pad=0.05, label="Z500 [m]", shrink=0.8)

# u as contour
axes[0, 0].contour(
    u_first.lon,
    u_first.lat,
    u_first,
    lw = 0.5,
    colors="black",
    levels=[level for level in u_levels_seq if level != 0],
    transform=ccrs.PlateCarree(),
)
axes[0, 1].contour(
    u_last.lon,
    u_last.lat,
    u_last,
    lw = 0.5,
    colors="black",
    levels=[level for level in u_levels_seq if level != 0],
    transform=ccrs.PlateCarree(),
)
axes[0, 2].contour(
    u_diff.lon,
    u_diff.lat,
    u_diff,
    lw = 0.5,
    colors="black",
    levels=[level for level in u_levels_div if level != 0],
    transform=ccrs.PlateCarree(),
)
# Second row for zg_std
# Plot zg_first_std
cf = axes[1, 0].contourf(
    zg_first_std.lon, zg_first_std.lat, zg_first_std, 
    cmap="viridis", 
    levels=zg_levels_seq,
    extend='both',
    transform=ccrs.PlateCarree()
)
axes[1, 0].set_title("first10 Z500 std", fontsize=15)
fig.colorbar(cf, ax=axes[1, 0], orientation='vertical', pad=0.05, label="Z500 [m]", shrink=0.8)

# Plot zg_last_std
cf = axes[1, 1].contourf(
    zg_last_std.lon, zg_last_std.lat, zg_last_std, 
    levels=zg_levels_seq,
    cmap="viridis", 
    extend='both',
    transform=ccrs.PlateCarree()
)
axes[1, 1].set_title("last10 Z500 std", fontsize=15)
fig.colorbar(cf, ax=axes[1, 1], orientation='vertical', pad=0.05, label="Z500 [m]", shrink=0.8)

# Plot zg_diff_std
cf = axes[1, 2].contourf(
    zg_diff_std.lon, zg_diff_std.lat, zg_diff_std, 
    levels=zg_levels_div,
    cmap="coolwarm", 
    extend='both',
    transform=ccrs.PlateCarree()
)
axes[1, 2].set_title("diff Z500 std", fontsize=15)
fig.colorbar(cf, ax=axes[1, 2], orientation='vertical', pad=0.05, label="Z500 [m]", shrink=0.8)

# u as contour
axes[1, 0].contour(
    u_first_nozon.lon,
    u_first_nozon.lat,
    u_first_nozon,
    lw = 0.5,
    colors="black",
    levels=[level for level in u_levels_seq_nonzon if level != 0],
    transform=ccrs.PlateCarree(),
)
axes[1, 1].contour(
    u_last_nozon.lon,
    u_last_nozon.lat,
    u_last_nozon,
    lw = 0.5,
    colors="black",
    levels=[level for level in u_levels_seq_nonzon if level != 0],
    transform=ccrs.PlateCarree(),
)
axes[1, 2].contour(
    u_diff_nozon.lon,
    u_diff_nozon.lat,
    u_diff_nozon,
    lw = 0.5,
    colors="black",
    levels=[level for level in u_levels_div if level != 0],
    transform=ccrs.PlateCarree(),
)

# third row for eof_rm
cf = axes[2, 0].contourf(
    eof_first_rm.lon,
    eof_first_rm.lat,
    eof_first_rm,
    cmap="RdBu_r",
    levels=zg_levels_div,
    extend='both',
    transform=ccrs.PlateCarree(),
)
axes[2,0].set_title("first10 NAO (18%)", fontsize=15)
fig.colorbar(cf, ax=axes[2, 0], orientation='vertical', pad=0.05, label="EOF_NAO [m]", shrink=0.8)

cf = axes[2, 1].contourf(
    eof_last_rm.lon,
    eof_last_rm.lat,
    eof_last_rm,
    cmap="RdBu_r",
    levels=zg_levels_div,
    extend='both',
    transform=ccrs.PlateCarree(),
)
axes[2,1].set_title("last10 NAO (24%)", fontsize=15)
fig.colorbar(cf, ax=axes[2, 1], orientation='vertical', pad=0.05, label="EOF_NAO[m]", shrink=0.8)

cf = axes[2, 2].contourf(
    eof_diff_rm.lon,
    eof_diff_rm.lat,
    eof_diff_rm,
    cmap="coolwarm",
    levels=zg_levels_div,
    extend='both',
    transform=ccrs.PlateCarree(),
)
axes[2,2].set_title("diff NAO", fontsize=15)
fig.colorbar(cf, ax=axes[2, 2], orientation='vertical', pad=0.05, label="EOF_NAO[m]", shrink=0.8)


for ax in axes.flat:
    ax.coastlines(lw = 0.5, color = 'grey7')
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

plt.tight_layout()
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/Z500_mean_std_NAO_nonzon.pdf")
# %%
