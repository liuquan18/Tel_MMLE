# %%
import xarray as xr
import numpy as np
import src.plots.utils as utils

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as mcolors
import proplot as pplt
from matplotlib.lines import Line2D

from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches


import src.plots.composite_plot as composite_plot
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import src.obs.era5_extreme_change as era5_extreme_change
import src.extreme.extreme_count_troposphere as ext_profile

#%%
def read_composite(
    model, var_name, fixed_pattern="decade_mpi", reduction="mean_same_number"
):
    """read composite data"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/"
    comp_name = f"plev_50000_{fixed_pattern}_first_JJA_JJA_first_last_{var_name}_composite_{reduction}.nc"
    composite = xr.open_dataset(f"{odir}{model}/composite/{comp_name}")
    if var_name == "ts":
        try:
            composite = composite.tsurf
        except AttributeError:
            composite = composite.ts
    elif var_name == "pr":
        try:
            composite = composite.pr
        except AttributeError:
            composite = composite.precip
    return composite


# read composite
def read_composite_rean(model, var_name, reduction="mean", group_size=40):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/composite/"
    composite_path = odir + "composite_mean_ts_40_withboot.nc"
    composite = xr.open_dataset(composite_path)
    return composite.__xarray_dataarray_variable__

# %%
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]

COMPOSITEs = {
    model: read_composite(model, var_name="pr", reduction='mean') for model in models
}

#%%
CR20_composite = read_composite_rean("CR20_allens", "pr")
#%%
# add to the dictionary
COMPOSITEs["20CR"] = CR20_composite


#%%
# change units from kg m-2 s-1 to mm/day
for model in COMPOSITEs.keys():
    COMPOSITEs[model] = COMPOSITEs[model] * 86400
#%%
models_plot = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "20CR"]
models_legend = [
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    "20CR (80)",
]



prec_cmap_seq = np.loadtxt(
    "/work/mh0033/m300883/High_frequecy_flow/data/colormaps-master/continuous_colormaps_rgb_0-1/prec_seq.txt"
)
prec_cmap_seq = mcolors.ListedColormap(prec_cmap_seq, name="prec_div")

prec_cmap_div = np.loadtxt(
    "/work/mh0033/m300883/High_frequecy_flow/data/colormaps-master/continuous_colormaps_rgb_0-1/prec_div.txt"
)
prec_cmap_div = mcolors.ListedColormap(prec_cmap_div, name="prec_div")
#%%
# %%

def plot_composite_single_ext(COMPOSITEs, models, axes, extr_type="pos", **kwargs):
    levels = kwargs.get("levels", np.arange(-1.2, 1.3, 0.3))
    for i, model in enumerate(models):  # cols for different models
        first = COMPOSITEs[model].sel(mode="NAO", period="first", extr_type=extr_type)
        last = COMPOSITEs[model].sel(mode="NAO", period="last", extr_type=extr_type)
        diff = COMPOSITEs[model].sel(mode="NAO", period="diff", extr_type=extr_type)
        diff_sig = COMPOSITEs[model].sel(
            mode="NAO", period="diff_sig", extr_type=extr_type
        )
        try:
            first = utils.erase_white_line(first)
            last = utils.erase_white_line(last)
            diff = utils.erase_white_line(diff)
            diff_sig = utils.erase_white_line(diff_sig)
        except ValueError:
            pass

        data_all = [
            first,
            last,
            diff,
        ]

        maps = []
        for j, data in enumerate(data_all):  # row for different data
            map = axes[j, i].contourf(
                data,
                x="lon",
                y="lat",
                levels=levels,
                extend="both",
                transform=ccrs.PlateCarree(),
                cmap = prec_cmap_div,
            )
            maps.append(map)
            axes[j, i].grid(color="grey7", linewidth=0.5)

        # significant area as hatches.
        if i < len(models) - 1:
            axes[2, i].contourf(
                diff_sig,
                levels=[-0.5, 0.5, 1.5],
                colors=["none", "none"],
                hatches=["", "xxxxx"],
            )

    return axes, maps

# %%

fig3, axes = pplt.subplots(
    space=0,
    width=180 / 25.4,
    wspace=0.2,
    hspace=0.2,
    proj="ortho",
    proj_kw=({"lon_0": -20, "lat_0": 60}),
    nrows=3,
    ncols=5,
)
axes.format(
    latlines=20,
    lonlines=30,
    color="grey7",
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
    leftlabels=["first", "last", "last - first"],
    toplabels=models_legend,
    toplabels_kw={"fontsize": 7, },
    leftlabels_kw={"fontsize": 7,},
)


axes, maps = plot_composite_single_ext(COMPOSITEs, models_plot, axes)
fig3.colorbar(
    maps[0],
    loc="b",
    pad=1,
    title="precipitation (mm/day)",
    width=0.1,
    shrink=1,
)

plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/precipitation_composite_pos.png")

# %%
fig4, axes = pplt.subplots(
    space=0,
    width=180 / 25.4,
    wspace=0.2,
    hspace=0.2,
    proj="ortho",
    proj_kw=({"lon_0": -20, "lat_0": 60}),
    nrows=3,
    ncols=5,
)
axes.format(
    latlines=20,
    lonlines=30,
    color="grey7",
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
    leftlabels=["first", "last", "last - first"],
    toplabels=models_legend,
    toplabels_kw={"fontsize": 7, },
    leftlabels_kw={"fontsize": 7,},
)


axes, maps = plot_composite_single_ext(COMPOSITEs, models, axes, 'neg')
fig4.colorbar(
    maps[0],
    loc="b",
    pad=1,
    title="precipitation (mm/day)",
    width=0.1,
    shrink=1,
)

plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/precipitation_composite_neg.png")
# %%
