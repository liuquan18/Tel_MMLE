# %%
import xarray as xr
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D
from src.mechnisms.mechisms import (
    read_jetStream,
    Jet_location,
    read_greenland_blocking,
    read_NAO_extremes,
)

import matplotlib.pyplot as plt

# %%

NAO_pos, NAO_neg = read_NAO_extremes("MPI_GE")
# %%
jet_stream = read_jetStream("MPI_GE")
jet_stream.load()
jet_loc = Jet_location(jet_stream)
jet_loc_north = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/mechnisms/jet_loc_north.nc"
).jet_loc
jet_loc_south = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/mechnisms/jet_loc_south.nc"
).jet_loc
# %%
blocking = read_greenland_blocking("MPI_GE")
blocking.load()
GB_pos = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/mechnisms/GB_pos.nc"
).GB
GB_neg = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/mechnisms/GB_neg.nc"
).GB
# %%
jet_clim = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/mechnisms/jet_loc_clim.nc"
).lat
jet_std = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/mechnisms/jet_loc_std.nc"
).lat

jet_upper = jet_clim + jet_std
jet_lower = jet_clim - jet_std
# %%
jet_clim_df = jet_clim.to_dataframe("jet_clim").reset_index()[["time", "jet_clim"]]
jet_upper_df = jet_upper.to_dataframe("jet_upper").reset_index()[["time", "jet_upper"]]
jet_lower_df = jet_lower.to_dataframe("jet_lower").reset_index()[["time", "jet_lower"]]


# %%
# select jet_loc based on NAO, nan to nan
jet_loc_NAO_pos = jet_loc.where(NAO_pos.notnull())
jet_loc_NAO_neg = jet_loc.where(NAO_neg.notnull())
# %%
jet_loc_north_NAO_pos = jet_loc_north.where(NAO_pos.notnull())
jet_loc_north_NAO_neg = jet_loc_north.where(NAO_neg.notnull())

jet_loc_south_NAO_pos = jet_loc_south.where(NAO_pos.notnull())
jet_loc_south_NAO_neg = jet_loc_south.where(NAO_neg.notnull())
# %%
GB_clim = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/mechnisms/GB_clim.nc"
).var156
GB_std = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/mechnisms/GB_std.nc"
).var156

GB_upper = GB_clim + GB_std
GB_lower = GB_clim - GB_std

# %%
GB_clim_df = GB_clim.to_dataframe("GB_clim").reset_index()[["time", "GB_clim"]]
GB_upper_df = GB_upper.to_dataframe("GB_upper").reset_index()[["time", "GB_upper"]]
GB_lower_df = GB_lower.to_dataframe("GB_lower").reset_index()[["time", "GB_lower"]]

# %%
blocking_NAO_pos = blocking.where(NAO_pos.notnull())
blocking_NAO_neg = blocking.where(NAO_neg.notnull())

GB_pos_NAO_pos = GB_pos.where(NAO_pos.notnull())
GB_pos_NAO_neg = GB_pos.where(NAO_neg.notnull())

GB_neg_NAO_pos = GB_neg.where(NAO_pos.notnull())
GB_neg_NAO_neg = GB_neg.where(NAO_neg.notnull())

# %%
# decadal count
NAO_pos_dec = (
    NAO_pos.resample(time="10Y", closed="left").count(dim=("time")).sum(dim="ens")
)
NAO_neg_dec = (
    NAO_neg.resample(time="10Y", closed="left").count(dim=("time")).sum(dim="ens")
)

# mean for absolute value (composite mean)
jet_loc_NAO_pos_dec = (
    jet_loc_NAO_pos.resample(time="10Y", closed="left").mean(dim="time").mean(dim="ens")
)
jet_loc_NAO_neg_dec = (
    jet_loc_NAO_neg.resample(time="10Y", closed="left")
    .mean(dim=("time"))
    .mean(dim="ens")
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

jet_loc_south_NAO_pos_dec = (
    jet_loc_south_NAO_pos.resample(time="10Y", closed="left")
    .count(dim=("time"))
    .sum(dim="ens")
)
jet_loc_south_NAO_neg_dec = (
    jet_loc_south_NAO_neg.resample(time="10Y", closed="left")
    .count(dim=("time"))
    .sum(dim="ens")
)

blocking_NAO_pos_dec = (
    blocking_NAO_pos.resample(time="10Y", closed="left")
    .mean(dim=("time"))
    .mean(dim="ens")
)
blocking_NAO_neg_dec = (
    blocking_NAO_neg.resample(time="10Y", closed="left")
    .mean(dim=("time"))
    .mean(dim="ens")
)

GB_pos_NAO_pos_dec = (
    GB_pos_NAO_pos.resample(time="10Y", closed="left")
    .count(dim=("time"))
    .sum(dim="ens")
)
GB_pos_NAO_neg_dec = (
    GB_pos_NAO_neg.resample(time="10Y", closed="left")
    .count(dim=("time"))
    .sum(dim="ens")
)

GB_neg_NAO_pos_dec = (
    GB_neg_NAO_pos.resample(time="10Y", closed="left")
    .count(dim=("time"))
    .sum(dim="ens")
)
GB_neg_NAO_neg_dec = (
    GB_neg_NAO_neg.resample(time="10Y", closed="left")
    .count(dim=("time"))
    .sum(dim="ens")
)

# %%
# to dataframe and concat along the column
NAO_pos_df = NAO_pos_dec.to_dataframe("NAO_pos").reset_index()[["time", "NAO_pos"]]
NAO_neg_df = NAO_neg_dec.to_dataframe("NAO_neg").reset_index()[["time", "NAO_neg"]]
jet_loc_NAO_pos_df = jet_loc_NAO_pos_dec.to_dataframe("jet_loc_compmean").reset_index()[
    ["time", "jet_loc_compmean"]
]
jet_loc_NAO_neg_df = jet_loc_NAO_neg_dec.to_dataframe("jet_loc_compmean").reset_index()[
    ["time", "jet_loc_compmean"]
]

jet_loc_north_NAO_pos_df = jet_loc_north_NAO_pos_dec.to_dataframe(
    "jet_north"
).reset_index()[["time", "jet_north"]]
jet_loc_north_NAO_neg_df = jet_loc_north_NAO_neg_dec.to_dataframe(
    "jet_north"
).reset_index()[["time", "jet_north"]]


jet_loc_south_NAO_pos_df = jet_loc_south_NAO_pos_dec.to_dataframe(
    "jet_south"
).reset_index()[["time", "jet_south"]]
jet_loc_south_NAO_neg_df = jet_loc_south_NAO_neg_dec.to_dataframe(
    "jet_south"
).reset_index()[["time", "jet_south"]]

blocking_NAO_pos_df = blocking_NAO_pos_dec.to_dataframe("GB_compmean").reset_index()[
    ["time", "GB_compmean"]
]
blocking_NAO_neg_df = blocking_NAO_neg_dec.to_dataframe("GB_compmean").reset_index()[
    ["time", "GB_compmean"]
]

GB_pos_NAO_pos_df = GB_pos_NAO_pos_dec.to_dataframe("GB_above").reset_index()[
    ["time", "GB_above"]
]
GB_pos_NAO_neg_df = GB_pos_NAO_neg_dec.to_dataframe("GB_above").reset_index()[
    ["time", "GB_above"]
]

GB_neg_NAO_pos_df = GB_neg_NAO_pos_dec.to_dataframe("GB_below").reset_index()[
    ["time", "GB_below"]
]
GB_neg_NAO_neg_df = GB_neg_NAO_neg_dec.to_dataframe("GB_below").reset_index()[
    ["time", "GB_below"]
]


# %%
# concatenate for positive NAO
NAO_pos_all_df = (
    NAO_pos_df.join(jet_loc_NAO_pos_df.set_index("time"), on="time")
    .join(jet_loc_north_NAO_pos_df.set_index("time"), on="time")
    .join(jet_loc_south_NAO_pos_df.set_index("time"), on="time")
    .join(blocking_NAO_pos_df.set_index("time"), on="time")
    .join(GB_pos_NAO_pos_df.set_index("time"), on="time")
    .join(GB_neg_NAO_pos_df.set_index("time"), on="time")
    .join(jet_clim_df.set_index("time"), on="time")
    .join(jet_upper_df.set_index("time"), on="time")
    .join(jet_lower_df.set_index("time"), on="time")
    .join(GB_clim_df.set_index("time"), on="time")
    .join(GB_upper_df.set_index("time"), on="time")
    .join(GB_lower_df.set_index("time"), on="time")
)

# %%
NAO_neg_all_df = (
    NAO_neg_df.join(jet_loc_NAO_neg_df.set_index("time"), on="time")
    .join(jet_loc_north_NAO_neg_df.set_index("time"), on="time")
    .join(jet_loc_south_NAO_neg_df.set_index("time"), on="time")
    .join(blocking_NAO_neg_df.set_index("time"), on="time")
    .join(GB_pos_NAO_neg_df.set_index("time"), on="time")
    .join(GB_neg_NAO_neg_df.set_index("time"), on="time")
    .join(jet_clim_df.set_index("time"), on="time")
    .join(jet_upper_df.set_index("time"), on="time")
    .join(jet_lower_df.set_index("time"), on="time")
    .join(GB_clim_df.set_index("time"), on="time")
    .join(GB_upper_df.set_index("time"), on="time")
    .join(GB_lower_df.set_index("time"), on="time")
)

# %%
NAO_pos_all_df["NAO_phase"] = "positive"
NAO_neg_all_df["NAO_phase"] = "negative"
# %%
# create column called decade, infer from time
NAO_pos_all_df["decade"] = NAO_pos_all_df["time"].dt.year // 10 * 10
NAO_neg_all_df["decade"] = NAO_neg_all_df["time"].dt.year // 10 * 10


# %%
# correspondence between NAO and background state
fig, axes = plt.subplots()

axes.plot(
    NAO_pos_all_df["decade"],
    (NAO_pos_all_df["jet_north"] / NAO_pos_all_df["NAO_pos"]),
    color='C1', linestyle='-'
)
axes.plot(
    NAO_neg_all_df["decade"],
    (NAO_neg_all_df["jet_north"] / NAO_neg_all_df["NAO_neg"]),
    color='C0', linestyle='-'
)
axes.plot(
    NAO_pos_all_df["decade"],
    (NAO_pos_all_df["GB_above"] / NAO_pos_all_df["NAO_pos"]),
    color='C1', linestyle='--'
)
axes.plot(
    NAO_neg_all_df["decade"],
    (NAO_neg_all_df["GB_above"] / NAO_neg_all_df["NAO_neg"]),
    color='C0', linestyle='--'
)
# Create custom legend

# Custom legend elements
legend_elements = [
    Line2D([0], [0], color='C1', linestyle='-', label='Positive NAO'),
    Line2D([0], [0], color='C1', linestyle='-', label='Negative NAO'),
    Line2D([0], [0], color='k',linestyle='-', label='eddy-driven jet stream'),
    Line2D([0], [0],color='k', linestyle='--', label='greenland blocking'),
]

# Add the legend to the plot
axes.legend(handles=legend_elements, loc='upper left')

axes.set_title('Correspondence between NAO and background state')

plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/NAO_background_state.png", dpi = 300)

# %%
NAO_all_df = pd.concat([NAO_pos_all_df, NAO_neg_all_df], ignore_index=True)

# %%
# new column 'extreme_count', value equals to 'NAO_pos' if NAO_phase is positive, else 'NAO_neg'
NAO_all_df["extreme_count"] = NAO_all_df.apply(
    lambda x: x["NAO_pos"] if x["NAO_phase"] == "positive" else x["NAO_neg"], axis=1
)
# %%
# new column, called jet_north/south, value equals to 'jet_north' if NAO_phase is positive, else -1 * 'jet_south'
NAO_all_df["jet_north_south"] = NAO_all_df.apply(
    lambda x: x["jet_north"] if x["NAO_phase"] == "positive" else -1 * x["jet_south"],
    axis=1,
)
# %%
fig, axes = plt.subplots(2, 1, figsize=(8, 10))
ax1 = axes[0]
scatter1 = sns.scatterplot(
    data=NAO_all_df,
    x="jet_north",
    y="NAO_pos",
    size="decade",
    sizes=(10, 500),
    legend="brief",
    ax=ax1,
    alpha=0.7,
    color = 'C1'
)
scatter2 = sns.scatterplot(
    data=NAO_all_df,
    x="jet_north",
    y="NAO_neg",
    size="decade",
    sizes=(10, 500),
    legend=False,
    ax=ax1,
    alpha=0.7,
    color = 'C0'
)
ax2 = axes[1]
sns.scatterplot(
    data=NAO_all_df,
    x="GB_above",
    y="NAO_pos",
    size="decade",
    sizes=(10, 500),
    legend="brief",
    ax=ax2,
    alpha = 0.7,
    color = 'C1'
)
sns.scatterplot(
    data=NAO_all_df,
    x="GB_above",
    y="NAO_neg",
    size="decade",
    sizes=(10, 500),
    legend=False,
    ax=ax2,
    alpha = 0.7,
    color = 'C0'
)

ax1.set_xlabel("count of jet stream location (> 1.5 std)")
ax1.set_ylabel("count of NAO extremes (> 1.5 std)")

ax2.set_xlabel("count of blocking (> 1.5 std)")
ax2.set_ylabel("count of NAO extremes (> 1.5 std)")

plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/NAO_GB_JET_count.png", dpi = 300)
# %%
fig, axes = plt.subplots(2, 1, figsize=(8, 10))
ax1 = axes[0]
scatter1 = sns.scatterplot(
    data=NAO_all_df,
    x="jet_loc_compmean",
    y="NAO_pos",
    hue="decade",
    size="decade",
    sizes=(200, 200),
    legend=False,
    ax=ax1,
    palette="flare",
)
scatter2 = sns.scatterplot(
    data=NAO_all_df,
    x="jet_loc_compmean",
    y="NAO_neg",
    hue="decade",
    size="decade",
    sizes=(200, 200),
    legend=False,
    ax=ax1,
    marker="x",
    palette="flare",
)
ax2 = axes[1]
sns.scatterplot(
    data=NAO_all_df,
    x="GB_compmean",
    y="NAO_pos",
    hue="decade",
    size="decade",
    sizes=(200, 200),
    legend="brief",
    ax=ax2,
    palette="flare",
)
sns.scatterplot(
    data=NAO_all_df,
    x="GB_compmean",
    y="NAO_neg",
    hue="decade",
    size="decade",
    sizes=(200, 200),
    legend=False,
    ax=ax2,
    marker="x",
    palette="flare",
)
# Move the legend to the bottom right
ax2.legend(loc="lower right")

ax1.set_xlabel(r"Eddy-driven jet stream location ($\degree$N)")
ax1.set_ylabel("NAO extremes (> 1.5 std)")

ax2.set_xlabel("Greenland Blocking Index proxy (km)")
ax2.set_ylabel("NAO extremes (> 1.5 std)")

plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/NAO_GB_JET_compmean.png", dpi = 300)


# %%


# %%
fig, axes = plt.subplots(2,1, figsize=(10, 10))

ax1 = axes[0]
sns.scatterplot(
    data=NAO_all_df,
    x="decade",
    y='jet_loc_compmean',
    hue="NAO_phase",
    legend="brief",
    palette="flare",
    size="extreme_count",
    sizes=(10, 500),
    ax=ax1,
)

sns.lineplot(data=NAO_all_df, x="decade", y="jet_clim", color="gray", ax=ax1)

sns.lineplot(data=NAO_all_df, x="decade", y="jet_upper", color="gray", ax=ax1)

sns.lineplot(data=NAO_all_df, x="decade", y="jet_lower", color="gray", ax=ax1)

ax2 = axes[1]
sns.scatterplot(
    data=NAO_all_df,
    x="decade",
    y="GB_compmean",
    hue="NAO_phase",
    legend="brief",
    palette="flare",
    size="extreme_count",
    sizes=(10, 500),
    ax=ax2,
)

sns.lineplot(data=NAO_all_df, x="decade", y="GB_clim", color="gray", ax=ax2)

sns.lineplot(data=NAO_all_df, x="decade", y="GB_upper", color="gray", ax=ax2)

sns.lineplot(data=NAO_all_df, x="decade", y="GB_lower", color="gray", ax=ax2)
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/temporal_NAO_JET_GB.png", dpi = 300)
#%%

fig, axes = plt.subplots(2,1, figsize=(10, 10))

ax1 = axes[0]
sns.scatterplot(
    data=NAO_all_df,
    x="extreme_count",
    y='jet_north',
    hue="decade",
    style="NAO_phase",
    legend="brief",
    palette="flare",
    size="extreme_count",
    sizes=(10, 500),
    ax=ax1,
)

ax2 = axes[1]
sns.scatterplot(
    data=NAO_all_df,
    x="extreme_count",
    y="GB_above",
    hue="decade",
    style="NAO_phase",
    legend="brief",
    palette="flare",
    size="extreme_count",
    sizes=(10, 500),
    ax=ax2,
)

# move legend to the outside
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# 45 degree line
ax1.plot([150, 250], [150, 250], ls="--", c=".3")
ax2.plot([150, 250], [150, 250], ls="--", c=".3")

plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/NAO_GB_JET_extreme_count_vs.png", dpi = 300)

# %%
sns.scatterplot(
    data=NAO_all_df,
    x="jet_loc_compmean",
    y="GB_compmean",
    hue="NAO_phase",
    legend="brief",
    palette="flare",
    size="extreme_count",
    sizes=(20, 200),
)
# %%
sns.scatterplot(
    data=NAO_all_df,
    x="jet_north_south",
    y="GB_above",
    hue="NAO_phase",
    legend="brief",
    palette="flare",
    size="extreme_count",
    sizes=(20, 200),
)
# %%
sns.scatterplot(
    data=NAO_all_df,
    x="jet_north_south",
    y="GB_above",
    style="NAO_phase",
    legend=False,
    palette="flare",
    size="extreme_count",
    hue="decade",
    sizes=(20, 500),
    alpha=0.7,
)


# %%
jet_loc_df = jet_loc.to_dataframe("jet_loc").reset_index()[["ens", "time", "jet_loc"]]
bloocking_df = blocking.to_dataframe("GB").reset_index()[["ens", "time", "GB"]]
jet_blocking_df = jet_loc_df.join(
    bloocking_df.set_index(["ens", "time"]), on=["ens", "time"]
)
# %%
jet_loc_extreme_NAO = jet_loc.where((NAO_pos.notnull()) | (NAO_neg.notnull()))
blocking_extreme_NAO = blocking.where((NAO_pos.notnull()) | (NAO_neg.notnull()))

jet_loc_extreme_NAO_df = jet_loc_extreme_NAO.to_dataframe("jet_loc").reset_index()[
    ["ens", "time", "jet_loc"]
]
blocking_extreme_NAO_df = blocking_extreme_NAO.to_dataframe("GB").reset_index()[
    ["ens", "time", "GB"]
]

jet_blocking_extreme_NAO_df = jet_loc_extreme_NAO_df.join(
    blocking_extreme_NAO_df.set_index(["ens", "time"]), on=["ens", "time"]
)

# %%
jet_loc_non_extreme_NAO = jet_loc.where((NAO_pos.isnull()) & (NAO_neg.isnull()))
blocking_non_extreme_NAO = blocking.where((NAO_pos.isnull()) & (NAO_neg.isnull()))

jet_loc_non_extreme_NAO_df = jet_loc_non_extreme_NAO.to_dataframe(
    "jet_loc"
).reset_index()[["ens", "time", "jet_loc"]]
blocking_non_extreme_NAO_df = blocking_non_extreme_NAO.to_dataframe("GB").reset_index()[
    ["ens", "time", "GB"]
]

jet_blocking_non_extreme_NAO_df = jet_loc_non_extreme_NAO_df.join(
    blocking_non_extreme_NAO_df.set_index(["ens", "time"]), on=["ens", "time"]
)


# %%
fig, ax = plt.subplots()
sns.scatterplot(
    data=NAO_all_df[NAO_all_df["decade"].isin(np.arange(1850, 1910, 10))],
    x="jet_loc_compmean",
    y="GB_compmean",
    style="NAO_phase",
    legend=False,
    palette="flare",
    size="extreme_count",
    hue="decade",
    sizes=(20, 500),
    alpha=0.7,
    ax=ax,
    edgecolor="black",  # Add edge color to the markers
    linewidth=1,  # Set the width of the edge
)

sns.scatterplot(
    data=NAO_all_df,
    x="jet_loc_compmean",
    y="GB_compmean",
    style="NAO_phase",
    legend="brief",
    palette="flare",
    size="extreme_count",
    hue="decade",
    sizes=(20, 500),
    alpha=0.7,
    ax=ax,
)


sns.lineplot(
    x="jet_clim", y="GB_clim", data=NAO_all_df, color="k", label="climatology", ax=ax
)


ax.set_xlabel(r"Eddy-driven jet stream location ($\degree$N)")
ax.set_ylabel("Greenland Blocking Index proxy (km)")

ax.set_xlim(40.5, 60.5)

# move the legend to the outside using bbox_to_anchor
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/NAO_GB_JET.png")


# %%
fig, ax = plt.subplots()
sns.scatterplot(
    data=NAO_all_df[NAO_all_df["decade"].isin(np.arange(1850, 1910, 10))],
    x="jet_loc_compmean",
    y="GB_compmean",
    style="NAO_phase",
    legend=False,
    palette="flare",
    size="extreme_count",
    hue="decade",
    sizes=(20, 500),
    alpha=0.7,
    ax=ax,
    edgecolor="black",  # Add edge color to the markers
    linewidth=1,  # Set the width of the edge
)

# Apply a linearly fitted line to the scatter plot
sns.regplot(
    data=NAO_all_df[NAO_all_df["NAO_phase"] == "negative"],
    x="jet_loc_compmean",
    y="GB_compmean",
    scatter=False,
    ax=ax,
    color="blue",
    line_kws={"label": "Linear Fit"},
    ci=None,  # Disable uncertainty shading
)

sns.scatterplot(
    data=NAO_all_df,
    x="jet_loc_compmean",
    y="GB_compmean",
    style="NAO_phase",
    legend="brief",
    palette="flare",
    size="extreme_count",
    hue="decade",
    sizes=(20, 500),
    alpha=0.7,
    ax=ax,
)

sns.regplot(
    data=NAO_all_df[NAO_all_df["NAO_phase"] == "positive"],
    x="jet_loc_compmean",
    y="GB_compmean",
    scatter=False,
    ax=ax,
    color="blue",
    line_kws={"label": "Linear Fit"},
    ci=None,  # Disable uncertainty shading
)


global_warming = sns.lineplot(
    x="jet_clim", y="GB_clim", data=NAO_all_df, color="k", label="climatology", ax=ax
)

non_extreme_first = sns.lineplot(
    x="jet_loc",
    y="GB",
    data=jet_blocking_non_extreme_NAO_df[
        jet_blocking_non_extreme_NAO_df["time"].dt.year.isin(np.arange(1850, 1910))
    ],
    color="k",
    ax=ax,
    alpha=0.2,
)

non_extreme_last = sns.lineplot(
    x="jet_loc",
    y="GB",
    data=jet_blocking_non_extreme_NAO_df[
        jet_blocking_non_extreme_NAO_df["time"].dt.year.isin(np.arange(2040, 2100))
    ],
    color="r",
    ax=ax,
    alpha=0.2,
)
# extremes only
extreme_first = sns.lineplot(
    x="jet_loc",
    y="GB",
    data=jet_blocking_extreme_NAO_df[
        jet_blocking_extreme_NAO_df["time"].dt.year.isin(np.arange(1850, 1980))
    ],
    color="k",
    ax=ax,
    alpha=0.2,
)

extreme_last = sns.lineplot(
    x="jet_loc",
    y="GB",
    data=jet_blocking_extreme_NAO_df[
        jet_blocking_extreme_NAO_df["time"].dt.year.isin(np.arange(2070, 2100))
    ],
    color="r",
    ax=ax,
    alpha=0.2,
)


ax.set_xlabel(r"Eddy-driven jet stream location ($\degree$N)")
ax.set_ylabel("Greenland Blocking Index proxy (km)")

ax.set_xlim(40.5, 60.5)

# move the legend to the outside using bbox_to_anchor
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
# plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/NAO_GB_JET.png")

# %%
