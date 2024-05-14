# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import src.plots.composite_plot as composite_plot
import proplot as pplt
import matplotlib.patches as patches
from matplotlib.widgets import PolygonSelector
# %%
ex = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/composite/plev_50000_decade_mpi_first_JJA_JJA_first_last_ts_composite_mean_same_number.nc"
)
# %%
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "20CR"]
models_legend = [
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    "20CR(80)",
]

fig3, axes = pplt.subplots(
    space=0,
    width=180 / 25.4,
    wspace=0.2,
    hspace=0.2,
    proj="ortho",
    proj_kw=({"lon_0": -20, "lat_0": 60}),
    nrows=3,
    ncols=1,
)
axes.format(
    latlines=20,
    lonlines=30,
    color="grey7",
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
    leftlabels=["first", "last", "last - first"],
    toplabels_kw={"fontsize": 7},
    leftlabels_kw={"fontsize": 7},
)

comps = {"MPI_GE": ex.tsurf}
axes, maps = composite_plot.plot_composite_single_ext(
    comps, models=["MPI_GE"], axes=axes
)
ax = axes[1, 0]
# y = [(397938.71034780145, -950631.363608636), (1326462.3678259999, -1746508.7843042351), (3360371.33182586, 1171708.4249129593), (2343416.849825926, 1967585.8456085622)]
# y.append(y[0]) #repeat the first point to create a 'closed loop'
y = [
[-19.8693420853922,53.53381413361307],
[8.431439164607797,41.104236521168126],
[29.297654068837293,49.9264460488725],
[63.099407914607795,60.75382198335459],
[26.712689164607795,72.49843690158606],
[-4.791851741217467,65.96430130824967],
[-19.8693420853922,53.53381413361307]
]

p = patches.Polygon(y, edgecolor = 'k', facecolor = 'none',transform=ccrs.PlateCarree())
ax.add_patch(p)

#%%
# create a boolean mask of the grid points that are inside the polygon
mask = ex.lon.apply(lambda x: p.contains_point(x[::-1]), axis=1)

# get the indices of the grid points that are inside the polygon
indices = np.where(mask)

# select the grid points that are inside the polygon
lonlat_inside_polygon = ex.lonlat.isel(lon=indices[0], lat=indices[1])

# print the longitude and latitude that is covered by the polygon
print(lonlat_inside_polygon)
#%%
plt.ion()    # interactive mode is on
fig3.show()
plt.show()

def onselect(data_input):
    print(data_input)
PS = PolygonSelector(axes[1,0], onselect)
a = input()    # prevent window from closing when execution is done

print("Click on the figure to create a polygon.")
print("Press the 'esc' key to start a new polygon.")
print("Try holding the 'shift' key to move all of the vertices.")
print("Try holding the 'ctrl' key to move a single vertex.")

print(PS.verts)
# %%
