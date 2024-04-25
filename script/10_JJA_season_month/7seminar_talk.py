#%%
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import proplot as pplt
import matplotlib.pyplot as plt
import matplotlib as mpl
#%%
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
plt.rcParams["hatch.linewidth"] = 0.3
plt.rcParams["text.color"] = "white"
plt.rcParams["axes.labelcolor"] = "white"
plt.rcParams["xtick.color"] = "white"
plt.rcParams["ytick.color"] = "white"
plt.rcParams["axes.edgecolor"] = "white"
plt.rcParams["figure.facecolor"] = "black"

# %%
eof_result = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CR20/EOF_result/first_plev550_eof.nc")
# %%
pattern = eof_result.eof.sel(mode = 'NAO').isel(decade = 0)

#%%
pc = eof_result.pc.sel(mode = 'NAO')
#%%
# reconstruct feild from pattern and pc
reconstruct = pattern * pc
# %%
fig1 = pplt.figure(figsize=(60 / 25.4, 60 / 25.4), sharex=False, sharey=False)
fig1.format(
    abc=False,
    facecolor="black",
)
# gs = pplt.GridSpec(
#     ncols=3,
#     nrows=1,
#     wspace=2,
#     wratios=(1, 1, 1),
# )
## Spatial pattern and index distribution
ax1 = fig1.add_subplot( proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax1.set_facecolor("white")

ax1.contourf(
    reconstruct.lon,
    reconstruct.lat,
    pattern,
    levels=np.arange(-30, 31, 5),
    extend = 'both',
    transform=ccrs.PlateCarree(),
    cmap = "RdBu_r",
)
ax1.format(
    lonlines=20,
    latlines=30,
    coast=True,
    coastlinewidth=0.5,
    coastcolor="charcoal",
)

plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/workshop/positive_NAO_pattern.pdf", facecolor=fig1.get_facecolor())
# %%
import scipy.stats as stats

# generate random data with normal distribution
data = np.random.normal(0, 1, 1000)

# convert data to percentiles
percentiles = stats.rankdata(data, "average") / len(data)

fig2 = pplt.figure(figsize=(60 / 25.4, 60 / 25.4), sharex=False, sharey=False)


# plot the histogram just lines
ax2 = fig2.add_subplot()
ax2.hist(percentiles, bins=30, histtype='step', color='white')

# plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/workshop/normal_distribution.pdf", facecolor=fig2.get_facecolor())
# %%
