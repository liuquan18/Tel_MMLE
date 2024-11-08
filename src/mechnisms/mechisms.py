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
def read_greenland_blocking(model):
    GreenlandBlocking = []
    gb_dir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/GB_"
    for month in ["Jun", "Jul", "Aug"]:
        all_ens_lists = sorted(glob.glob(gb_dir + month + "/*.nc"))
        gb = xr.open_mfdataset(all_ens_lists, combine="nested", concat_dim="ens")
        GreenlandBlocking.append(gb)

    GreenlandBlocking = xr.concat(GreenlandBlocking, dim="time")
    # sort by time
    GreenlandBlocking = GreenlandBlocking.sortby("time")
    # exclude lon dim
    GreenlandBlocking = GreenlandBlocking.var156.isel(lat = 0, lon = 0, plev = 0)
    # exclude data of 2100
    GreenlandBlocking = GreenlandBlocking.sel(time=slice("1850", "2099"))
    # chagne from m to km
    GreenlandBlocking = GreenlandBlocking / 1000
    return GreenlandBlocking

