# %%
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import os

from src.mechnisms.mechisms import *

import pandas as pd
import seaborn as sns
import numpy as np

#%%
def decade_df(NAO_dec, jet_dec, GB_dec, phase):
    df = pd.DataFrame(
        {
            "decade": NAO_dec.time.dt.year// 10 * 10,
            "extreme_count": NAO_dec.values,
            "jet_loc_north": jet_dec.values,
            "GB_above": GB_dec.values,
            "phase": phase,
        }
    )
    return df

# %%
model = "MK36"
# %%
def save_results(model):
    jet_loc, jet_loc_clim, jet_loc_std, jet_loc_north, jet_north_decade = process_jet(model)
    GB, GB_clim, GB_std,GB_above, GB_above_decade = process_GB(model)

    NAO_pos, NAO_neg = read_NAO_extremes(model)
    
    jet_loc_north_NAO_pos, jet_loc_north_NAO_neg = NAO_correspond(NAO_pos, NAO_neg, jet_loc_north)
    GB_above_NAO_pos, GB_above_NAO_neg = NAO_correspond(NAO_pos, NAO_neg, GB_above)

    
    # decadal count
    NAO_pos_dec = (
        NAO_pos.resample(time="10Y", closed="left").count(dim=("time")).sum(dim="ens")
    )
    NAO_neg_dec = (
        NAO_neg.resample(time="10Y", closed="left").count(dim=("time")).sum(dim="ens")
    )

    
    jet_loc_north_NAO_pos_dec = (
        jet_loc_north_NAO_pos.resample(time="10Y", closed="left")
        .count(dim=("time"))
        .sum(dim="ens")
    )
    jet_loc_north_NAO_neg_dec = (
        jet_loc_north_NAO_neg.resample(time="10Y", closed="left")
        .count(dim=("time"))
        .sum(dim="ens")
    )

    
    GB_above_NAO_pos_dec = (
        GB_above_NAO_pos.resample(time="10Y", closed="left")
        .count(dim=("time"))
        .sum(dim="ens")
    )

    GB_above_NAO_neg_dec = (
        GB_above_NAO_neg.resample(time="10Y", closed="left")
        .count(dim=("time"))
        .sum(dim="ens")
    )

    NAO_pos_df = decade_df(NAO_pos_dec, jet_loc_north_NAO_pos_dec, GB_above_NAO_pos_dec, "positive")
    NAO_neg_df = decade_df(NAO_neg_dec, jet_loc_north_NAO_neg_dec, GB_above_NAO_neg_dec, "negative")

    # save results
    to_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/"
    if not os.path.exists(to_dir):
        os.makedirs(to_dir)

    NAO_pos_df.to_csv(to_dir + "NAO_pos_df.csv")
    NAO_neg_df.to_csv(to_dir + "NAO_neg_df.csv")
# %%
for model in ["MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]:
    print(model)
    save_results(model)
# %%
save_results('CanESM2')
# %%
