# %%
import xarray as xr
import pandas as pd
import seaborn as sns
from src.mechnisms.mechisms import (
    read_jetStream,
    Jet_location,
    read_greenland_blocking,
    read_NAO_extremes,
    to_dataframe,
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
# select jet_loc based on NAO, nan to nan
jet_loc_NAO_pos = jet_loc.where(NAO_pos.notnull())
jet_loc_NAO_neg = jet_loc.where(NAO_neg.notnull())
# %%
jet_loc_north_NAO_pos = jet_loc_north.where(NAO_pos.notnull())
jet_loc_north_NAO_neg = jet_loc_north.where(NAO_neg.notnull())

jet_loc_south_NAO_pos = jet_loc_south.where(NAO_pos.notnull())
jet_loc_south_NAO_neg = jet_loc_south.where(NAO_neg.notnull())
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
)

# %%
NAO_neg_all_df = (
    NAO_neg_df.join(jet_loc_NAO_neg_df.set_index("time"), on="time")
    .join(jet_loc_north_NAO_neg_df.set_index("time"), on="time")
    .join(jet_loc_south_NAO_neg_df.set_index("time"), on="time")
    .join(blocking_NAO_neg_df.set_index("time"), on="time")
    .join(GB_pos_NAO_neg_df.set_index("time"), on="time")
    .join(GB_neg_NAO_neg_df.set_index("time"), on="time")
)

#%%
NAO_pos_all_df['NAO_phase'] = 'positive'
NAO_neg_all_df['NAO_phase'] = 'negative'
#%%
# create column called decade, infer from time
NAO_pos_all_df['decade'] = NAO_pos_all_df['time'].dt.year // 10 * 10
NAO_neg_all_df['decade'] = NAO_neg_all_df['time'].dt.year // 10 * 10
#%%
NAO_all_df = pd.concat([NAO_pos_all_df, NAO_neg_all_df], ignore_index=True)

# %%
fig, ax = plt.subplots()
scatter1 = sns.scatterplot(data=NAO_all_df, x="jet_north", y='NAO_pos', hue='decade', size='decade', sizes=(20, 200), legend="full", ax=ax)
scatter2 = sns.scatterplot(data=NAO_all_df, x="jet_north", y='NAO_neg', hue='decade', size='decade', sizes=(20, 200), legend="full", ax=ax, marker='x')

# add colorbar at the bottom for hue
norm = plt.Normalize(NAO_all_df['decade'].min(), NAO_all_df['decade'].max())
sm = plt.cm.ScalarMappable(cmap="deep", norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', pad=0.2)
cbar.set_label('Decade')


# %%
fig, ax = plt.subplots()
sns.scatterplot(data=NAO_all_df, x="GB_above", y = 'NAO_pos', hue = 'decade', size = 'decade',sizes = (20, 200), legend="full", ax = ax)
sns.scatterplot(data=NAO_all_df, x="GB_above", y = 'NAO_neg', hue = 'decade', size = 'decade',sizes = (20, 200), legend="full", ax = ax, marker = 'x')
# %%
sns.scatterplot(data=NAO_all_df, x="jet_south", y = 'NAO_neg', hue = 'decade', size = 'decade',sizes = (20, 200), legend="full", marker = 'x')

# %%
