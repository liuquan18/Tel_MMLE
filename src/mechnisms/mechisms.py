#%%
import xarray as xr
import glob

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

