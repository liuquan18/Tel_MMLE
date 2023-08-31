#%%
import pandas as pd
import matplotlib as mpl
import proplot as pplt
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#%%
import src.obs.era5_extreme_change as era5_extreme_change
#%%


mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams['font.family'] = 'sans-serif'
plt.rcParams['hatch.linewidth'] = 0.3

# black background
mpl.rcParams['axes.facecolor'] = 'black'
mpl.rcParams['figure.facecolor'] = 'black'
mpl.rcParams['savefig.facecolor'] = 'black'

# white font color
mpl.rcParams['text.color'] = 'white'
mpl.rcParams['axes.labelcolor'] = 'white'
mpl.rcParams['axes.edgecolor'] = 'white'
mpl.rcParams['xtick.color'] = 'white'
mpl.rcParams['ytick.color'] = 'white'
plt.rcParams['xtick.color'] = 'white'
plt.rcParams['ytick.color'] = 'white'


#######################
# prepare daily data
#%%
url = 'https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.nao.cdas.z500.19500101_current.csv'
df = pd.read_csv(url)
# %%

# Combine year, month, and day columns into a single datetime column
df['time'] = pd.date_range(start='1950-01-01', end='2023-08-29', freq='D')

# Drop the original year, month, and day columns
df.drop(['year', 'month', 'day'], axis=1, inplace=True)

# Rename the remaining column to 'nao_index'
df.rename(columns={'nao_index_cdas': 'nao_index'}, inplace=True)
# %%
# set time as index
df = df.set_index('time')
# %%
# plot the line plot of df
df_jja_2023 = df.loc['2023-06':'2023-08']

#####################
# plot
# %%
fig, ax = plt.subplots(figsize=(6, 4))
df_jja_2023.plot(ax = ax,kind='line',color='w',linewidth=2)

ax.set_ylim(-3.2,3.2)
ax.set_yticks([-3, -1.5, 0, 1.5, 3])

ax.set_xlim('2023-06-01','2023-08-31')

# hline at 0
ax.axhline(y=0, color='w', linestyle='dotted',linewidth=1)

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
# save
plt.savefig('/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/nao_index_2023.pdf',bbox_inches='tight',dpi=300)

#%%
#####################
# load month data
ERA5 = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/ERA5/EOF_result/plev_50000_1940_2022_ERA5_all.nc")
ERA5 = ERA5.pc.sel(mode = 'NAO').squeeze()

ERA5_nodec = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/ERA5/EOF_result/plev_50000_1940_2022_ERA5_no_dec_all.nc")
ERA5_nodec = ERA5_nodec.pc.sel(mode = 'NAO').squeeze()

# %%
fig, ax = plt.subplots(figsize=(6, 4))
df_jja_2023.plot(ax = ax,kind='line',color='w',linewidth=2)

ax.set_ylim(-3.2,3.2)
ax.set_yticks([-3, -1.5, 0, 1.5, 3])

ax = era5_extreme_change.plot_era_nao_index(
    ERA5,
    ERA5_nodec,
    ax,
)
# %%
