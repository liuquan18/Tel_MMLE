# %%
import pandas as pd
import matplotlib as mpl
import proplot as pplt
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import matplotlib.ticker as ticker

# %%
import src.obs.era5_extreme_change as era5_extreme_change

#%%
import importlib
importlib.reload(era5_extreme_change)
# %%


mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
plt.rcParams["hatch.linewidth"] = 0.3

# black background
mpl.rcParams["axes.facecolor"] = "black"
mpl.rcParams["figure.facecolor"] = "black"
mpl.rcParams["savefig.facecolor"] = "black"

# white font color
mpl.rcParams["text.color"] = "white"
mpl.rcParams["axes.labelcolor"] = "white"
mpl.rcParams["axes.edgecolor"] = "white"
mpl.rcParams["xtick.color"] = "white"
mpl.rcParams["ytick.color"] = "white"
plt.rcParams["xtick.color"] = "white"
plt.rcParams["ytick.color"] = "white"


#######################
# prepare daily data
# %%
url = "https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.nao.cdas.z500.19500101_current.csv"
df = pd.read_csv(url)
# %%

# Combine year, month, and day columns into a single datetime column
df["time"] = pd.date_range(start="1950-01-01", end="2023-08-29", freq="D")

# Drop the original year, month, and day columns
df.drop(["year", "month", "day"], axis=1, inplace=True)

# Rename the remaining column to 'nao_index'
df.rename(columns={"nao_index_cdas": "nao_index"}, inplace=True)
# %%
# set time as index
df = df.set_index("time")
# %%
# plot the line plot of df
df_jja_2023 = df.loc["2023-06":"2023-08"]

#####################
# plot
# %%
fig, ax = plt.subplots(figsize=(6, 4))
df_jja_2023.plot(ax=ax, kind="line", color="w", linewidth=2)

ax.set_ylim(-3.2, 3.2)
ax.set_yticks([-3, -1.5, 0, 1.5, 3])

ax.set_xlim("2023-06-01", "2023-08-31")

# hline at 0
ax.axhline(y=0, color="w", linestyle="dotted", linewidth=1)

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
# change the ticks
ax.tick_params(
    axis="x",
    which="major",
    direction="out",
    pad=2,
    labelsize=7,
    labelcolor="white",
)

ax.tick_params(
    axis="y",
    which="major",
    direction="out",
    pad=2,
    labelsize=7,
    labelcolor="white",
)

# set null to minor ticks
ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())

# set xmajor ticks as the first, tenth, twentieth, and thirtieth day of the month
ax.xaxis.set_major_locator(mpl.dates.DayLocator(bymonthday=[1, 10, 20]))
ax.xaxis.set_major_formatter(mpl.dates.DateFormatter("%d %b"))

# no legend and title
ax.legend().set_visible(False)
ax.set_ylabel("NAO standard deviation")
ax.set_xlabel("")
# save
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/nao_index_2023.pdf",
    bbox_inches="tight",
    dpi=300,
)
#%%
df_jja_2023 = df.loc["2023-06":"2023-08"]

# monthly mean of df_jja_2023
df_jja_2023_monthly = df_jja_2023.resample("M").mean()
#%%
df_jja_2023_monthly.index = pd.to_datetime(['2023-06-01', '2023-07-01', '2023-08-01'])
df_jja_2023_monthly = df_jja_2023_monthly.to_xarray()
df_jja_2023_monthly = df_jja_2023_monthly.rename({'index':'time'}).nao_index

# %%
#####################
# load month data
ERA5 = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/ERA5/EOF_result/plev_50000_1940_2022_ERA5_all.nc"
)
ERA5 = ERA5.pc.sel(mode="NAO").drop_vars(('mode','plev')).squeeze()

# add also the data for 2023 from the daily data mean
ERA5 = xr.concat([ERA5, df_jja_2023_monthly],dim = 'time')


ERA5_nodec = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/ERA5/EOF_result/plev_50000_1940_2022_ERA5_no_dec_all.nc"
)
ERA5_nodec = ERA5_nodec.pc.sel(mode="NAO").drop_vars('mode').squeeze()

# %%
fig, ax = plt.subplots(figsize=(6, 4))

ax.set_ylim(-3.2, 3.2)
ax.set_yticks([-3, -1.5, 0, 1.5, 3])

nao = ERA5
NEW_nao = ERA5_nodec.sel(time=slice("1941", "2022"))


# hline at y= 1.5 and -1.5
threshod = 1.5
ax.axhline(y=threshod, color="w", linestyle="dotted", linewidth=1)
ax.axhline(y=-1 * threshod, color="w", linestyle="dotted", linewidth=1)
ax.axhline(y=0, color="w", linestyle="dotted", linewidth=1)

first_pos_org, first_neg_org, last_pos_org, last_neg_org = era5_extreme_change.count_extreme(
    nao, threshod=threshod
)
first_pos_new, first_neg_new, last_pos_new, last_neg_new = era5_extreme_change.count_extreme(
    NEW_nao, threshod=threshod
)

ax.plot(nao.values, color="white", alpha=1, lw=1.5, label="NAO")
# ax.plot(NEW_nao.values, color="red", lw=0.5, label="NAO no decade")

# vline at x = 1981
xmin, xmax = ax.get_xlim()
xmid = (xmin + xmax) / 2
# ax.axvline(x=xmid, color="g", linestyle="--")
ax.set_yticks([-3, -1.5, 0, 1.5, 3])


ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
# change the ticks
ax.tick_params(
    axis="x",
    which="major",
    direction="out",
    pad=2,
    labelsize=7,
    labelcolor="white",
)

ax.tick_params(
    axis="y",
    which="major",
    direction="out",
    pad=2,
    labelsize=7,
    labelcolor="white",
)

# set null to minor ticks
ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())

# # put the count as text on the plot
# # pos
# ax.text(
#     0.20,
#     0.99,
#     f"1941-1981: {first_pos_org}",
#     transform=ax.transAxes,
#     fontsize=7,
#     color="gray",
#     alpha=0.8,
#     verticalalignment="top",
# )
# ax.text(
#     0.55,
#     0.99,
#     f"1982-2022: {last_pos_org}",
#     transform=ax.transAxes,
#     fontsize=7,
#     color="gray",
#     alpha=0.8,
#     verticalalignment="top",
# )

# ax.text(
#     0.20,
#     0.93,
#     f"1941-1981: {first_pos_new}",
#     transform=ax.transAxes,
#     fontsize=7,
#     color = 'red',
#     alpha=0.8,
#     verticalalignment="top",
# )
# ax.text(
#     0.55,
#     0.93,
#     f"1982-2022: {last_pos_new}",
#     transform=ax.transAxes,
#     fontsize=7,
#     color = 'red',
#     alpha=0.8,
#     verticalalignment="top",
# )

# # neg
# ax.text(
#     0.20,
#     0.15,
#     f"1941-1981: {first_neg_org}",
#     transform=ax.transAxes,
#     fontsize=7,
#     color="gray",
#     alpha=0.8,
#     verticalalignment="top",
# )
# ax.text(
#     0.55,
#     0.15,
#     f"1982-2022: {last_neg_org}",
#     transform=ax.transAxes,
#     fontsize=7,
#     color="gray",
#     alpha=0.8,
#     verticalalignment="top",
# )

# ax.text(
#     0.20,
#     0.10,
#     f"1941-1981: {first_neg_new}",
#     transform=ax.transAxes,
#     fontsize=7,
#     color = 'red',
#     alpha=0.8,
#     verticalalignment="top",
# )
# ax.text(
#     0.55,
#     0.10,
#     f"1982-2022: {last_neg_new}",
#     transform=ax.transAxes,
#     fontsize=7,
#     color = 'red',
#     alpha=0.8,
#     verticalalignment="top",
# )
# ax.set_title("")
# ax.legend(loc="top", fontsize=7, frameon=False)
ax.xaxis.set_major_locator(
    ticker.FixedLocator(np.arange(-1, 246, 30))
)  # every 10 years (JJA)
ax.xaxis.set_major_formatter(plt.FuncFormatter(era5_extreme_change.format_year_summer))
# ax.set_xlim(-1, 246)
ax.set_ylabel("NAO standard deviation")
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/nao_index_1941_2022_threshod.pdf",
)
# %%
