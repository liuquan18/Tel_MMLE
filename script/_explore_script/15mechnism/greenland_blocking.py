# %%
import xarray as xr
import pandas as pd
import seaborn as sns
import glob
import os
import matplotlib.pyplot as plt

from src.mechnisms.mechisms import read_greenland_blocking
# %%
# %%
def GreenlandBlocking_climatology(GreenlandBlocking, stat='mean'):
    if stat == 'mean':
        GreenlandBlocking_clim = GreenlandBlocking.resample(time="10Y", closed='left').mean(dim=('time', 'ens'))
    else:
        GreenlandBlocking_clim = GreenlandBlocking.resample(time="10Y", closed='left').std(dim=('time', 'ens'))

    return GreenlandBlocking_clim

# %%
def decade_sub_GreenlandBlocking(GreenlandBlocking, GreenlandBlocking_clim):

    # Initialize an empty array to store the results
    GreenlandBlocking_anomaly = GreenlandBlocking.copy()

    # Iterate over each decade
    for i in range(len(GreenlandBlocking_clim.time)):
        # Get the start and end of the decade
        end_year = GreenlandBlocking_clim.time[i].dt.year.item()
        start_year = end_year - 9

        # Select the years within the current decade
        decade_years = GreenlandBlocking.sel(time=slice(f"{start_year}", f"{end_year}"))

        # Subtract the climatology of the decade from each year in the decade
        GreenlandBlocking_anomaly.loc[dict(time=decade_years.time)] = decade_years - GreenlandBlocking_clim.isel(time=i)

    return GreenlandBlocking_anomaly

#%%
def decade_jet_NS(jet_loc, jet_loc_clim, jet_loc_std, scale = 1, fix_clim = False):

    jet_loc_north = jet_loc.copy()
    jet_loc_south = jet_loc.copy()

    for i in range (len(jet_loc_clim.time)):
        end_year = jet_loc_clim.time[i].dt.year.item()
        start_year = end_year - 9

        decade_years = jet_loc.sel(time=slice(f"{start_year}", f"{end_year}"))

        if fix_clim:
            clim_mean = jet_loc_clim.isel(time=0)
            clim_std = jet_loc_std.isel(time=0)
        else:
            clim_mean = jet_loc_clim.isel(time=i)
            clim_std = jet_loc_std.isel(time=i)

        jet_loc_north.loc[dict(time=decade_years.time)] = decade_years.where(decade_years > clim_mean + scale * clim_std)
        jet_loc_south.loc[dict(time=decade_years.time)] = decade_years.where(decade_years < clim_mean - scale * clim_std)

    return jet_loc_north, jet_loc_south

# %%
def blocking_process(model):
    GB = read_greenland_blocking(model)
    GB.load()

    GB.name = 'GB'

    GB_clim = GreenlandBlocking_climatology(GB)

    GB_std = GreenlandBlocking_climatology(GB, stat='std')
    GB = GB.drop_vars(('lat','lon', 'plev'))
    GB_clim = GB_clim.drop_vars(('lat','lon', 'plev'))
    GB_std = GB_std.drop_vars(('lat','lon', 'plev'))
    GB_pos, GB_neg = decade_jet_NS(GB, GB_clim, GB_std, scale = 1, fix_clim = True)



    GB_pos_decade = GB_pos.resample(time="10Y", closed='left').count(dim=('time')).sum(dim = 'ens')



    # save results
    clim_to_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/"

    # create dir if not exist
    if not os.path.exists(clim_to_dir):
        os.makedirs(clim_to_dir)

    GB_clim.to_netcdf(clim_to_dir + "GB_clim.nc")
    GB_std.to_netcdf(clim_to_dir + "GB_std.nc")

    GB_pos_decade.to_netcdf(clim_to_dir + "GB_pos_decade.nc")

    return GB_clim, GB_std, GB_pos_decade


# %%
GB_clim, GB_std, GB_pos_decade = blocking_process('CanESM2')
fig, ax = plt.subplots()
GB_clim.plot(ax=ax, label = 'GB_clim', color = 'black')
# fill between 1 std
ax.fill_between(GB_clim.time, GB_clim - GB_std, GB_clim + GB_std, alpha=0.3, color = 'gray', label = '1 std')

ax_twin = ax.twinx()
GB_pos_decade.plot(ax=ax_twin, label = 'count of GB above 1std', color = 'red')

# ax_twin.set_ylim(-1200, 2700)
ax.set_ylabel('GB proxy (km)')
ax.set_title('Greenland Blocking Index')
# put legend together
lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax_twin.get_legend_handles_labels()
ax.legend(lines + lines2, labels + labels2, loc='upper left', frameon=False)
# plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/GB_index.png")
# %%
for model in ['MPI_GE', "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]:
    blocking_process(model)
# %%
