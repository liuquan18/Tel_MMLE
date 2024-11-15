# %%
import xarray as xr
import pandas as pd
import seaborn as sns
import glob
import matplotlib.pyplot as plt
import os
from src.mechnisms.mechisms import read_jetStream, Jet_location
# %%

# %%
def Jet_climatology(JetStream, stat = 'mean'):
    if stat == 'mean':
        JetStream_clim = JetStream.resample(time="10Y", closed = 'left').mean(dim = ('time','ens'))
    else:
        JetStream_clim = JetStream.resample(time="10Y", closed = 'left').std(dim = ('time','ens'))

    return JetStream_clim


# %%
def decade_sub(JetStream, JetStream_clim):

    # Initialize an empty array to store the results
    JetStream_anomaly = JetStream.copy()

    # Iterate over each decade
    for i in range(len(JetStream_clim.time)):
        # Get the start and end of the decade
        end_year = JetStream_clim.time[i].dt.year.item()
        start_year = end_year - 9

        # Select the years within the current decade
        decade_years = JetStream.sel(time=slice(f"{start_year}", f"{end_year}"))

        # Subtract the climatology of the decade from each year in the decade
        JetStream_anomaly.loc[dict(time=decade_years.time)] = decade_years - JetStream_clim.isel(time=i)

    return JetStream_anomaly
#%%
def decade_jet_NS(jet_loc, jet_loc_clim, jet_loc_std, scale = 1.5, fix_clim = False):

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
def zonal_wind_process(model):
    jet_stream = read_jetStream(model)
    jet_stream = jet_stream.compute()

    jet_loc = Jet_location(jet_stream)

    jet_loc.name = 'jet_loc'

    jet_loc_clim = Jet_climatology(jet_loc)
    jet_loc_clim = jet_loc_clim.drop_vars('lon')

    jet_loc_std = Jet_climatology(jet_loc, stat = 'std')
    jet_loc_std = jet_loc_std.drop_vars('lon')



    jet_loc_north, jet_loc_south = decade_jet_NS(jet_loc, jet_loc_clim, jet_loc_std, fix_clim=True)


    jet_north_decade = jet_loc_north.resample(time="10Y", closed = 'left').count().sum(dim = 'ens')


    # save results
    clim_to_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/"

    # create dir if not exist
    if not os.path.exists(clim_to_dir):
        os.makedirs(clim_to_dir)

    jet_loc_clim.to_netcdf(clim_to_dir + 'jet_loc_clim.nc')
    jet_loc_std.to_netcdf(clim_to_dir + 'jet_loc_std.nc')

    jet_north_decade.to_netcdf(clim_to_dir + 'jet_north_decade.nc')

    return jet_loc_clim, jet_loc_std, jet_north_decade


# %%
jet_loc_clim, jet_loc_std, jet_north_decade = zonal_wind_process('CanESM2')


fig, ax = plt.subplots()
jet_loc_clim.plot(ax=ax, label='Climatology',x = 'time', add_legend = False)
# clean the tile of ax
ax.set_title('')

ax.fill_between(jet_loc_clim.time, jet_loc_clim - jet_loc_std, jet_loc_clim + jet_loc_std, color='gray', alpha=0.5, label = 'std')
ax.set_ylim(42,57)
ax_twin = ax.twinx()
jet_north_decade.drop_vars('lon').plot(ax = ax_twin, color = 'red', label = 'count of jet Norther than 1 std')
# ax_twin.set_ylim(-800, 1000)

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax_twin.get_legend_handles_labels()
ax.legend(lines + lines2, labels + labels2, loc='upper left', frameon=False)

ax.set_ylabel('Latitude')
ax_twin.set_ylabel('count')
ax.set_title('Eddy-driven jet stream location')
# plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/jet_stream.png", dpi = 300)

# %%
for model in ["MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]:
    zonal_wind_process(model)
# %%
