# %%
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.pyplot as plt
# %%
from src.mechnisms.mechisms import *


# %%
model = "CESM1_CAM5"

# %%
jet_stream = read_jetStream(model)
jet_stream = jet_stream.load()
# %%
jet_loc = Jet_location(jet_stream)
jet_loc_clim = decade_climatology(jet_loc)
jet_loc_clim = jet_loc_clim.drop_vars("lon")
# %%
jet_loc_std = decade_climatology(jet_loc, stat='std')
jet_loc_std = jet_loc_std.drop_vars("lon")
# %%
jet_loc_north, jet_loc_south = decade_classify(jet_loc, jet_loc_clim, jet_loc_std, scale = 1.5, fix_clim=True)

jet_loc_north.drop_vars("lon")
jet_loc_south.drop_vars("lon")
# %%
jet_north_decade = jet_loc_north.resample(time="10Y", closed = 'left').count().sum(dim = 'ens')


# %%
GB = read_greenland_blocking(model)
GB.load()
# %%
GB_clim = decade_climatology(GB)
GB_clim = GB_clim.drop_vars(["lon", "lat", "plev"])
# %%
GB_std = decade_climatology(GB, stat='std')
GB_std = GB_std.drop_vars(["lon", "lat", "plev"])
# %%
GB_above, GB_below = decade_classify(GB, GB_clim, GB_std, scale = 1.5, fix_clim=True)
GB_above = GB_above.drop_vars(["lon", "lat", "plev"])
GB_below = GB_below.drop_vars(["lon", "lat", "plev"])
# %%
GB_above_decade = GB_above.resample(time="10Y", closed = 'left').count().sum(dim = 'ens')
GB_below_decade = GB_below.resample(time="10Y", closed = 'left').count().sum(dim = 'ens')
# %%
def decade_corr(x,y):
    x = x.stack(com = ('time','ens'))
    y = y.stack(com = ('time','ens'))
    corr = xr.corr(x, y, dim='com')
    return corr
# %%
jet_GB_corr = jet_loc.resample(time = '10Y',closed="left" ).apply(
    lambda x: decade_corr(x, GB.sel(time = x.time))
)
# %%
jet_GB_corr = jet_GB_corr.drop_vars(('lon', 'lat', 'plev'))

# %%
fig = plt.figure(figsize=(6, 8))
gs = fig.add_gridspec(3, 1)

# Left panel for time series
ax1 = fig.add_subplot(gs[0])

jet_loc_clim.plot(ax=ax1, label='Climatology', x='time', add_legend=False, color='black')
# clean the tile of ax
ax1.set_title('')

ax1.fill_between(jet_loc_clim.time, jet_loc_clim - jet_loc_std, jet_loc_clim + jet_loc_std, color='gray', alpha=0.5, label='std')
ax1.set_ylim(42, 57)
ax1.grid(False)  # Remove gridlines
ax_twin1 = ax1.twinx()
jet_north_decade.drop_vars('lon').plot(ax=ax_twin1, color='red', label='count of jet Norther than 1.5 std')
ax_twin1.grid(False)  # Remove gridlines

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax_twin1.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper left', frameon=False)

ax1.set_ylabel('Latitude')
ax_twin1.set_ylabel('count')
ax1.set_title('Eddy-driven jet stream location')

ax2 = fig.add_subplot(gs[1])
GB_clim.plot(ax=ax2, label='Climatology', x='time', add_legend=False, color='black')
# clean the tile of ax
ax2.set_title('')

ax2.fill_between(GB_clim.time, GB_clim - GB_std, GB_clim + GB_std, color='gray', alpha=0.5, label='std')
ax2.set_ylim(5.45, 5.7)
ax2.grid(False)  # Remove gridlines
ax_twin2 = ax2.twinx()
GB_above_decade.plot(ax=ax_twin2, color='red', label='count of GB above 1.5 std')
ax_twin2.grid(False)  # Remove gridlines

lines, labels = ax2.get_legend_handles_labels()
lines2, labels2 = ax_twin2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper left', frameon=False)

ax2.set_ylabel('GB proxy (km)')
ax_twin2.set_ylabel('count')
ax2.set_title('Greenland Blocking Index')

ax3 = fig.add_subplot(gs[2])
jet_GB_corr.plot(ax=ax3, color = 'k')
ax3.set_ylabel("Correlation")
ax3.set_title("Jet Stream location and Blocking correlation")
# without grids
ax3.grid(False)
ax3.set_ylim(-0.5, 0)
plt.tight_layout()

# plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/jet_GB_index_evolve.png")
# %%
