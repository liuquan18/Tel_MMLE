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
def corresponding_df(NAO_pos, NAO_neg, jet_loc, GB):
    NAO_pos_df = (
        NAO_pos.to_dataframe("NAO_pos")
        .reset_index()
        .dropna(subset=["NAO_pos"])[["ens", "time", "NAO_pos"]]
    )
    NAO_pos_df["NAO_phase"] = "positive"

    NAO_neg_df = (
        NAO_neg.to_dataframe("NAO_neg")
        .reset_index()
        .dropna(subset=["NAO_neg"])[["ens", "time", "NAO_neg"]]
    )
    NAO_neg_df["NAO_phase"] = "negative"

    jet_loc_NAO_pos, jet_loc_NAO_neg = NAO_correspond(NAO_pos, NAO_neg, jet_loc)
    GB_NAO_pos, GB_NAO_neg = NAO_correspond(NAO_pos, NAO_neg, GB)

    # Jet Stream
    jet_loc_NAO_pos_df = (
        jet_loc_NAO_pos.to_dataframe("jet_loc")
        .reset_index()
        .dropna(subset=["jet_loc"])[["ens", "time", "jet_loc"]]
    )

    jet_loc_NAO_neg_df = (
        jet_loc_NAO_neg.to_dataframe("jet_loc")
        .reset_index()
        .dropna(subset=["jet_loc"])[["ens", "time", "jet_loc"]]
    )

    # Greenland Blocking
    GB_NAO_pos_df = (
        GB_NAO_pos.to_dataframe("GB")
        .reset_index()
        .dropna(subset=["GB"])[["ens", "time", "GB"]]
    )

    GB_NAO_neg_df = (
        GB_NAO_neg.to_dataframe("GB")
        .reset_index()
        .dropna(subset=["GB"])[["ens", "time", "GB"]]
    )

    # join pos

    NAO_pos_jet_GB = NAO_pos_df.join(
        jet_loc_NAO_pos_df.set_index(["ens", "time"]), on=["ens", "time"]
    ).join(GB_NAO_pos_df.set_index(["ens", "time"]), on=["ens", "time"])

    NAO_neg_jet_GB = NAO_neg_df.join(
        jet_loc_NAO_neg_df.set_index(["ens", "time"]), on=["ens", "time"]
    ).join(GB_NAO_neg_df.set_index(["ens", "time"]), on=["ens", "time"])

    NAO_pos_jet_GB["decade"] = NAO_pos_jet_GB["time"].dt.year // 10 * 10
    NAO_neg_jet_GB["decade"] = NAO_neg_jet_GB["time"].dt.year // 10 * 10

    return NAO_pos_jet_GB, NAO_neg_jet_GB


# %%
def split_period(
    NAO_pos_jet_GB, NAO_neg_jet_GB, jet_loc_clim, GB_clim, start=1850, end=2090
):

    NAO_pos_jet_GB_plot = NAO_pos_jet_GB[NAO_pos_jet_GB["decade"].isin([start, end])]
    # add column called period, value equals to 'first10' if decade is 1850, else 'last10'
    NAO_pos_jet_GB_plot["period"] = NAO_pos_jet_GB_plot["decade"].apply(
        lambda x: "first10" if x == start else "last10"
    )

    NAO_neg_jet_GB_plot = NAO_neg_jet_GB[NAO_neg_jet_GB["decade"].isin([start, end])]
    # add column called period, value equals to 'first10' if decade is 1850, else 'last10'
    NAO_neg_jet_GB_plot["period"] = NAO_neg_jet_GB_plot["decade"].apply(
        lambda x: "first10" if x == start else "last10"
    )

    jet_clim_first10 = jet_loc_clim.sel(time=str(start + 9))
    jet_clim_last10 = jet_loc_clim.sel(time=str(end + 9))

    GB_clim_first10 = GB_clim.sel(time=str(start + 9))
    GB_clim_last10 = GB_clim.sel(time=str(end + 9))
    return (
        NAO_pos_jet_GB_plot,
        NAO_neg_jet_GB_plot,
        jet_clim_first10,
        jet_clim_last10,
        GB_clim_first10,
        GB_clim_last10,
    )

# %%
# %%
# read data
def save_results(model):
    jet_loc, jet_loc_clim, jet_loc_std, jet_loc_north, jet_north_decade = process_jet(model)
    GB, GB_clim, GB_std,GB_above, GB_above_decade = process_GB(model)

    NAO_pos, NAO_neg = read_NAO_extremes(model)

    jet_GB_corr = correlate_jet_GB(jet_loc, GB)

    NAO_pos_jet_GB, NAO_neg_jet_GB = corresponding_df(NAO_pos, NAO_neg, jet_loc, GB)

    (
        NAO_pos_jet_GB_plot,
        NAO_neg_jet_GB_plot,
        jet_clim_first10,
        jet_clim_last10,
        GB_clim_first10,
        GB_clim_last10,
    ) = split_period(NAO_pos_jet_GB, NAO_neg_jet_GB, jet_loc_clim, GB_clim, 1950, 2090)

    # save results
    to_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/"
    if not os.path.exists(to_dir):
        os.makedirs(to_dir)

    NAO_pos_jet_GB.to_csv(to_dir + "NAO_pos_jet_GB_correspond.csv")
    NAO_neg_jet_GB.to_csv(to_dir + "NAO_neg_jet_GB_correspond.csv")

#%%
for model in ["MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]:
    print(model)
    save_results(model)



# %%
fig, ax = plt.subplots(figsize=(8, 8))
sns.kdeplot(
    data=NAO_pos_jet_GB_plot,
    x="jet_loc",
    y="GB",
    ax=ax,
    hue="period",
    fill=True,
    alpha=0.5,
    common_norm=True,
    legend=False,
)

sns.kdeplot(
    data=NAO_neg_jet_GB_plot,
    x="jet_loc",
    y="GB",
    ax=ax,
    hue="period",
    fill=False,
    alpha=0.7,
    common_norm=True,
    legend=False,
)


ax.set_xlabel(r"Eddy-driven jet stream location ($\degree$N)")
ax.set_ylabel("Greenland Blocking Index proxy (km)")

ax.set_title("extreme NAO under different background states")


# vline for jet_clim_first10 and jet_clim_last10
ax.axvline(
    jet_clim_first10, linestyle="--", label="jet climatology (1850-1859)", color="C0"
)
ax.axvline(
    jet_clim_last10, linestyle="--", label="jet climatology (2090-2099)", color="C1"
)
# hline for GB_clim_first10 and GB_clim_last10
ax.axhline(
    GB_clim_first10, linestyle="--", label="GB climatology (1850-1859)", color="C0"
)
ax.axhline(
    GB_clim_last10, linestyle="--", label="GB climatology (2090-2099)", color="C1"
)

# Create custom legend
legend_elements = [
    Line2D([0], [0], color="C0", lw=2, label="first10 (1850-1859)", linestyle="-"),
    Line2D([0], [0], color="C1", lw=2, label="last10 (2090-2099)", linestyle="-"),
    Patch(facecolor="C0", edgecolor="k", label="NAO (positive)", alpha=0.5),
    Patch(facecolor="none", edgecolor="k", label="NAO (negative)", alpha=0.5),
    Line2D([0], [0], color="C0", lw=2, label="first10 climatology", linestyle="--"),
    Line2D([0], [0], color="C1", lw=2, label="last10 climatology", linestyle="--"),
]


ax.legend(handles=legend_elements, loc="upper right")
plt.tight_layout()
# %%
