# %%
import xarray as xr
import pandas as pd
import seaborn as sns
import glob
import matplotlib.pyplot as plt
# %%


############## read NAO ###################
def read_NAO_extremes(
    model, standard="first", season="JJA", vertical_eof="ind", plev=50000
):
    eof_dir = (
        f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"
        + f"troposphere_{vertical_eof}_decade_"
        + standard
        + "_"
        + season
        + "_eof_result.nc"
    )

    eof_result = xr.open_dataset(eof_dir)

    pc = eof_result.pc.sel(mode="NAO", plev=plev)

    # extremes
    pos = pc.where(pc > 1.5)
    pos_df = pos.to_dataframe().reset_index()
    pos_df = pos_df.dropna(subset=["pc"])

    neg = pc.where(pc < -1.5)
    neg_df = neg.to_dataframe().reset_index()
    neg_df = neg_df.dropna(subset=["pc"])

    return pos_df, neg_df


# %%
NAO_pos, NAO_neg = read_NAO_extremes("MPI_GE")


# %%
def read_jetStream(model):
    JetStream = []
    jet_dir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/NA_EJ_"
    for month in ["Jun", "Jul", "Aug"]:
        all_ens_lists = sorted(glob.glob(jet_dir + month + "/*.nc"))
        jet = xr.open_mfdataset(all_ens_lists, combine="nested", concat_dim="ens")
        JetStream.append(jet)

    JetStream = xr.concat(JetStream, dim="time")
    # sort by time
    JetStream = JetStream.sortby("time")
    # exclude lon dim
    JetStream = JetStream.isel(lon=0).var131
    # exclude data of 2100
    JetStream = JetStream.sel(time=slice("1850", "2099"))
    return JetStream


# %%
def Jet_location(JetStream):
    jet_loc = JetStream.lat[JetStream.argmax(dim="lat")]
    return jet_loc


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
def decade_jet_NS(jet_loc, jet_loc_clim, jet_loc_std):
    jet_loc_north = jet_loc.copy()
    jet_loc_south = jet_loc.copy()

    for i in range (len(jet_loc_clim.time)):
        end_year = jet_loc_clim.time[i].dt.year.item()
        start_year = end_year - 9

        decade_years = jet_loc.sel(time=slice(f"{start_year}", f"{end_year}"))

        jet_loc_north.loc[dict(time=decade_years.time)] = decade_years.where(decade_years > jet_loc_clim.isel(time=i) + jet_loc_std.isel(time=i))
        jet_loc_south.loc[dict(time=decade_years.time)] = decade_years.where(decade_years < jet_loc_clim.isel(time=i) - jet_loc_std.isel(time=i))

    return jet_loc_north, jet_loc_south

# %%
jet_stream = read_jetStream("MPI_GE")
jet_stream = jet_stream.compute()
# %%
jet_loc = Jet_location(jet_stream)
# %%
jet_loc_clim = Jet_climatology(jet_loc)
#%%
jet_loc_std = Jet_climatology(jet_loc, stat = 'std')
#%%
fig, ax = plt.subplots()
jet_loc_clim.plot(ax=ax)
# jet_stream.isel(ens = 0).plot(x = 'time', ax=ax)
# fill between the std
ax.fill_between(jet_loc_clim.time, jet_loc_clim - jet_loc_std, jet_loc_clim + jet_loc_std, color='gray', alpha=0.5)
jet_loc_north.plot.scatter(x = 'time', ax = ax)
jet_loc_south.plot.scatter(x = 'time', ax = ax)
# %%
jet_loc_north, jet_loc_south = decade_jet_NS(jet_loc, jet_loc_clim, jet_loc_std)
# %%
# above 1 std clim
# %%
jet_north_decade = jet_loc_north.resample(time="10Y", closed = 'left').count().sum(dim = 'ens')
# %%
jet_south_decade = jet_loc_south.resample(time="10Y", closed = 'left').count().sum(dim = 'ens')
# %%
