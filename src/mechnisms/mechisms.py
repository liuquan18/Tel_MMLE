#%%
import xarray as xr
import glob

# %%

############## read NAO ###################
def read_NAO_extremes(
    model, 
):
    eof_dir = (f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result.nc"
    )

    eof_result = xr.open_dataset(eof_dir)

    pc = eof_result.pc.sel(mode="NAO")

    # extremes
    pos = pc.where(pc > 1.5)

    neg = pc.where(pc < -1.5)

    return pos, neg

def to_dataframe(arr, name = 'pc'):
    df = arr.to_dataframe(name).reset_index()
    df = df.dropna(subset=[name])
    return df
# %%
def read_jetStream(model):
    JetStream = []
    jet_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/NA_EJ_"
    for month in ["Jun", "Jul", "Aug"]:
        all_ens_lists = sorted(glob.glob(jet_dir + month + "/*.nc"))
        jet = xr.open_mfdataset(all_ens_lists, combine="nested", concat_dim="ens")
        JetStream.append(jet)

    JetStream = xr.concat(JetStream, dim="time")
    # sort by time
    JetStream = JetStream.sortby("time")
    # exclude lon dim
    JetStream = JetStream.isel(lon=0)

    try:
        JetStream = JetStream.var131
    except AttributeError:
        JetStream = JetStream.ua

    # time format
    try:
        JetStream["time"] = JetStream.indexes["time"].to_datetimeindex()
    except:
        pass

    # exclude data of 2100
    JetStream = JetStream.sel(time=slice("1850", "2099"))
    return JetStream


# %%
def Jet_location(JetStream):
    jet_loc = JetStream.lat[JetStream.argmax(dim="lat")]
    return jet_loc

# %%
def read_greenland_blocking(model):
    GreenlandBlocking = []
    gb_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/GB_"
    for month in ["Jun", "Jul", "Aug"]:
        all_ens_lists = sorted(glob.glob(gb_dir + month + "/*.nc"))
        gb = xr.open_mfdataset(all_ens_lists, combine="nested", concat_dim="ens")
        GreenlandBlocking.append(gb)

    GreenlandBlocking = xr.concat(GreenlandBlocking, dim="time")
    # sort by time
    GreenlandBlocking = GreenlandBlocking.sortby("time")
    # exclude lon dim
    GreenlandBlocking = GreenlandBlocking.isel(lat = 0, lon = 0, plev = 0)
    try:
        GreenlandBlocking = GreenlandBlocking.var156
    except AttributeError:
        try:
            GreenlandBlocking = GreenlandBlocking.zg
        except AttributeError:
            GreenlandBlocking = GreenlandBlocking.gph

    try:
        GreenlandBlocking["time"] = GreenlandBlocking.indexes["time"].to_datetimeindex()
    except:
        pass
    # exclude data of 2100
    GreenlandBlocking = GreenlandBlocking.sel(time=slice("1850", "2099"))
    # chagne from m to km
    GreenlandBlocking = GreenlandBlocking / 1000
    return GreenlandBlocking

# %%
def decade_climatology(var, stat='mean'):
    if stat == 'mean':
        var_clim = var.resample(time="10Y", closed='left').mean(dim=('time', 'ens'))
    else:
        var_clim = var.resample(time="10Y", closed='left').std(dim=('time', 'ens'))

    return var_clim
# %%
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
# %%
def decade_classify(var, var_clim, var_std, scale = 1, fix_clim = True):
    
        var_north = var.copy()
        var_south = var.copy()
    
        for i in range (len(var_clim.time)):
            end_year = var_clim.time[i].dt.year.item()
            start_year = end_year - 9
    
            decade_years = var.sel(time=slice(f"{start_year}", f"{end_year}"))
    
            if fix_clim:
                clim_mean = var_clim.isel(time=0)
                clim_std = var_std.isel(time=0)
            else:
                clim_mean = var_clim.isel(time=i)
                clim_std = var_std.isel(time=i)
    
            var_north.loc[dict(time=decade_years.time)] = decade_years.where(decade_years > clim_mean + scale * clim_std)
            var_south.loc[dict(time=decade_years.time)] = decade_years.where(decade_years < clim_mean - scale * clim_std)
    
        return var_north, var_south