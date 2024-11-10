# %%
import xarray as xr
import pandas as pd
import seaborn as sns
import glob
import matplotlib.pyplot as plt
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
jet_stream = read_jetStream("MPI_GE")
jet_stream = jet_stream.compute()
# %%
jet_loc = Jet_location(jet_stream)
# %%
jet_loc_clim = Jet_climatology(jet_loc)
jet_loc_clim = jet_loc_clim.drop_vars('lon')
#%%
jet_loc_std = Jet_climatology(jet_loc, stat = 'std')
jet_loc_std = jet_loc_std.drop_vars('lon')

# %%
jet_loc_north, jet_loc_south = decade_jet_NS(jet_loc, jet_loc_clim, jet_loc_std, fix_clim=True)

#%%
jet_loc_north.name = 'jet_loc'
jet_loc_south.name = 'jet_loc'

jet_loc_north.to_netcdf('/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/mechnisms/jet_loc_north.nc')
jet_loc_south.to_netcdf('/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/mechnisms/jet_loc_south.nc')
# %%
# above 1 std clim
# %%
jet_north_decade = jet_loc_north.resample(time="10Y", closed = 'left').count().sum(dim = 'ens')
# %%
jet_south_decade = jet_loc_south.resample(time="10Y", closed = 'left').count().sum(dim = 'ens')

#%%
fig, ax = plt.subplots()
jet_loc_clim.plot(ax=ax)
# jet_stream.isel(ens = 0).plot(x = 'time', ax=ax)
# fill between the std
ax.fill_between(jet_loc_clim.time, jet_loc_clim - jet_loc_std, jet_loc_clim + jet_loc_std, color='gray', alpha=0.5)
jet_loc_north.plot.scatter(x = 'time', ax = ax)
jet_loc_south.plot.scatter(x = 'time', ax = ax)
# %%
jet_north_decade.plot()
# %%
jet_south_decade.plot()
#%%

# %%
fig, ax = plt.subplots()
jet_loc_clim.plot(ax=ax, label='Climatology',x = 'time', add_legend = False)
# clean the tile of ax
ax.set_title('')

ax.fill_between(jet_loc_clim.time, jet_loc_clim - jet_loc_std, jet_loc_clim + jet_loc_std, color='gray', alpha=0.5, label = 'std')
ax.set_ylim(42,57)
ax_twin = ax.twinx()
jet_north_decade.drop_vars('lon').plot(ax = ax_twin, color = 'red', label = 'count of jet Norther than 1 std')
ax_twin.set_ylim(-800, 1000)

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax_twin.get_legend_handles_labels()
ax.legend(lines + lines2, labels + labels2, loc='upper left', frameon=False)

ax.set_ylabel('Latitude')
ax_twin.set_ylabel('count')
ax.set_title('Eddy-driven jet stream location')
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/jet_stream.png", dpi = 300)

# %%
jet_loc_clim_perc = (jet_loc_clim - jet_loc_clim.isel(time=0))/jet_loc_clim.isel(time=0) * 100
jet_loc_upper = jet_loc_clim + jet_loc_std
jet_loc_lower = jet_loc_clim - jet_loc_std

jet_loc_upper_perc = (jet_loc_upper - jet_loc_upper.isel(time=0))/jet_loc_upper.isel(time=0) * 100
jet_loc_lower_perc = (jet_loc_lower - jet_loc_lower.isel(time=0))/jet_loc_lower.isel(time=0) * 100

#%%
fig, ax = plt.subplots()
jet_loc_clim_perc.plot(ax=ax, label='Climatology')
ax.fill_between(jet_loc_clim_perc.time, jet_loc_lower_perc, jet_loc_upper_perc, color='gray', alpha=0.5, label = 'std')

ax.set_ylabel('Latitude anomaly (%)')
ax.set_title('eddy-driven Jet stream location')

plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/jet_stream_anomaly.png", dpi = 300)

#%%