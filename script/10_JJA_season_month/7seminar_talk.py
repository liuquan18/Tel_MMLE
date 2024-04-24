#%%
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import proplot as pplt


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
fig1 = pplt.figure(figsize=(180 / 25.4, 60 / 25.4), sharex=False, sharey=False)
fig1.format(
    abc=False,
    facecolor="black",
)
gs = pplt.GridSpec(
    ncols=3,
    nrows=1,
    wspace=2,
    wratios=(1, 1, 1),
)
## Spatial pattern and index distribution
ax1 = fig1.add_subplot( gs[0],proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax1.set_facecolor("white")

ax1.contourf(
    reconstruct.lon,
    reconstruct.lat,
    reconstruct.isel(time = 0),
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



# %%
import matplotlib.animation as animation

# Create a function that generates each frame
def update(i):
    ax1.clear()
    ax1.set_facecolor("white")
    ax1.contourf(
        reconstruct.lon,
        reconstruct.lat,
        reconstruct.isel(time=i),
        levels=np.arange(-30, 31, 5),
        extend='both',
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
    )
    ax1.format(
        lonlines=20,
        latlines=30,
        coast=True,
        coastlinewidth=0.5,
        coastcolor="charcoal",
    )

# Create the animation
ani = animation.FuncAnimation(fig1, update, frames=120)

# Save the animation
ani.save('/work/mh0033/m300883/Tel_MMLE/docs/source/plots/workshop/animation.mp4', writer='ffmpeg')
# %%
