# %%
import xarray as xr
import pandas as pd
import seaborn as sns
import glob
# %%

############## read NAO ###################
def read_NAO_extremes( model, standard="first", season="JJA",vertical_eof = 'ind', plev = 50000):
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
    pos_df = pos_df.dropna(subset=['pc'])

    
    neg = pc.where(pc < -1.5)
    neg_df = neg.to_dataframe().reset_index()
    neg_df = neg_df.dropna(subset=['pc'])

    return pos_df,neg_df
# %%
NAO_pos, NAO_neg = read_NAO_extremes('MPI_GE')

# %%
def read_jetStream(model):
    JetStream  = []
    jet_dir = "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/NA_EJ_"
    for month in ['Jun', 'Jul', 'Aug']:
        all_ens_lists = sorted(glob.glob(jet_dir + month + "/*.nc"))
        jet = xr.open_mfdataset(all_ens_lists, combine= 'nested', concat_dim = 'ens')
        JetStream.append(jet)

    JetStream = xr.concat(JetStream, dim = 'time')
    # sort by time
    JetStream = JetStream.sortby('time')
    # exclude lon dim
    JetStream = JetStream.isel(lon = 0)
    return JetStream
# %%
JetStream = read_jetStream('MPI_GE')
# %%
