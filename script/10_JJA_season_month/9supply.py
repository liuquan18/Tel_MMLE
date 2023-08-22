# %%
import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
import proplot as pplt
import seaborn as sns
import cartopy.crs as ccrs
import matplotlib.ticker as ticker
from matplotlib import lines as mlines
import matplotlib.ticker as ticker



from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
import matplotlib.patches as mpatches


import src.plots.composite_plot as composite_plot
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import src.plots.utils as utils
import src.obs.era5_extreme_change as era5_extreme_change

# %%
import importlib

importlib.reload(extplt)
importlib.reload(stat_overview)
importlib.reload(composite_plot)
importlib.reload(era5_extreme_change)


# %%
def read_composite(model, var_name,reduction = 'mean'):
    """read composite data"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/"
    comp_name = f"plev_50000_decade_mpi_first_JJA_JJA_first_last_{var_name}_composite_{reduction}.nc"
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
    elif var_name == "zg":
        try:
            composite = composite.zg
        except AttributeError:
            composite = composite.var156
    return composite

def read_all_models(variable):
    """read all models data"""
    models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    if variable == "eof_decade":
        all_model_data = {model: read_eof_decade(model) for model in models}
    elif variable == "eof_all":
        all_model_data = {
            model: read_eof_all(model) for model in models[1:]
        }  # no onepct here
        all_model_data["ERA5"] = read_eof_all("ERA5")  # also add ERA5
        all_model_data["ERA5_no_dec"] = read_eof_all("ERA5_no_dec") # also add ERA5_no_dec
    elif variable == "extrc":
        all_model_data = {model: read_extrc(model) for model in models}
    elif variable == "composite":
        models = models[1:]
        all_model_data = {
            model: read_composite(model, var_name="zg") for model in models
        }
    return all_model_data
# %%
COMPOSITEs = read_all_models("composite")
# %%

# %%
# Fig 3, composite plot of ts for positve extremes
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
models_legend = [
"MPI-GE (100)",
"CanESM2 (50)",
"CESM1-CAM5 (40)",
"MK3.6 (30)",
"GFDL-CM3 (20)",
]

fig3, axes = pplt.subplots(
    space=0,
    width = 180/25.4,
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
    color = 'grey7',
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
    leftlabels=["first10", "last10", "last10 - first10"],
    toplabels=models_legend,
    toplabels_kw = {"fontsize": 7},
    leftlabels_kw = {"fontsize": 7},
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

axes,maps = composite_plot.plot_composite_single_ext(COMPOSITEs, models, axes,levels = np.arange(-60,61,10))
fig3.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K",width = 0.1,shrink=1)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/composite_zg_pos.png",
    dpi=300,
    bbox_inches="tight",
)


# %%
# Fig 4, composite plot of ts for neg extremes
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
models_legend = [
"MPI-GE (100)",
"CanESM2 (50)",
"CESM1-CAM5 (40)",
"MK3.6 (30)",
"GFDL-CM3 (20)",
]

fig4, axes = pplt.subplots(
    space=0,
    width = 180/25.4,
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
    color = 'grey7',
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
    leftlabels=["first10", "last10", "last10 - first10"],
    toplabels=models_legend,
    toplabels_kw = {"fontsize": 7},
    leftlabels_kw = {"fontsize": 7},
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

axes,maps = composite_plot.plot_composite_single_ext(COMPOSITEs, models, axes,extr_type='neg',levels = np.arange(-60,61,10))
fig4.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K",width = 0.1,shrink=1)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/composite_zg_neg.png",
    dpi=300,
    bbox_inches="tight",
)

# %%
