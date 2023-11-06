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
######################
# read data
######################

# SMILEs
def read_eof_decade(model):
    """read eofs that is decomposed by decade"""
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"
    filename = "plev_50000_decade_mpi_first_JJA_eof_result.nc"
    ds = xr.open_dataset(odir + filename)
    ds = ds.sel(mode="NAO")
    return ds


def read_extrc(model):
    """read extreme counts"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_count/"
    filename = "plev_50000_decade_mpi_first_JJA_extre_counts.nc"
    ds = xr.open_dataset(odir + filename).pc
    return ds


def read_composite(model, var_name,reduction = 'mean_same_number'):
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
    return composite


def read_all_models(variable):
    """read all models data"""
    models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    if variable == "eof_decade":
        all_model_data = {model: read_eof_decade(model) for model in models}
    elif variable == "extrc":
        all_model_data = {model: read_extrc(model) for model in models}
    elif variable == "composite":
        models = models[1:] # no MPI_GE_onepct
        all_model_data = {
            model: read_composite(model, var_name="ts") for model in models
        }
    return all_model_data

#%%
# 20CR
# eof
def read_eof_rean(model,group_size = 40):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    first_eof_path = odir + f"EOF_result/first_{str(group_size)}_eof_std.nc"
    last_eof_path = odir + f"EOF_result/last_{str(group_size)}_eof_std.nc"

    first_eof = xr.open_dataset(first_eof_path)
    last_eof = xr.open_dataset(last_eof_path)
    return first_eof, last_eof

# read extreme counts
def read_extrc_rean(model, group_size = 40):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    first_extc_path = odir + f'first_{str(group_size)}_extc.nc'
    last_extc_path = odir + f'last_{str(group_size)}_extc.nc'

    first_extc = xr.open_dataset(first_extc_path)
    last_extc = xr.open_dataset(last_extc_path)
    return first_extc, last_extc

# read composite
def read_composite_rean(model, var_name,reduction = 'mean',group_size = 40):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/composite/"
    composite_path = odir + f'composite_{reduction}_{var_name}_{group_size}.nc'
    composite = xr.open_dataset(composite_path)
    return composite

# %%
######################
# utils functions
######################

def split_first_last(eof_result):
    times = eof_result.time
    years = np.unique(times.dt.year)
    first_years = years[:10]
    last_years = years[-10:]

    eof_first = eof_result.isel(decade=0).sel(
        time=eof_result["time.year"].isin(first_years)
    )
    eof_last = eof_result.isel(decade=-1).sel(
        time=eof_result["time.year"].isin(last_years)
    )
    return eof_first, eof_last


# %%
######################
# matplotlib global parameters
######################

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams['font.family'] = 'sans-serif'
plt.rcParams['hatch.linewidth'] = 0.3
# %%
# SMILEs
EOFs_decade = read_all_models("eof_decade")
EXTRCs = read_all_models("extrc")
COMPOSITEs = read_all_models("composite")

#%%
# 20CR
CR20_first_eof, CR20_last_eof = read_eof_rean('CR20_allens')
CR20_first_extc, CR20_last_extc = read_extrc_rean('CR20_allens')
CR20_composite = read_composite_rean('CR20_allens','ts')


# %%
# Fig 1, spatial pattern chagne, index distribution change, and extreme lines
# set the fig size as 180mm x 150mm
# set the font size as 7pt
fig1 = pplt.figure(figsize=(180 / 25.4, 180 / 25.4),sharex=False,sharey=False)
fig1.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
    
)
models_legend = [
    "MPI_GE_onepct (100)",
    "MPI-GE_Hist+RCP8.5 (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    ]

gs = pplt.GridSpec(
    ncols=4,
    nrows=2,
    # wspace=(5, 0.5),
    # hspace=0.5,
    hratios=[1, 1],
    wratios=[1,1, 1,1],
)

# MPI_GE
ax1 = fig1.add_subplot(gs[0, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax2 = fig1.add_subplot(gs[0,1])

# 20CR
ax3 = fig1.add_subplot(gs[1, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax4 = fig1.add_subplot(gs[1,1])

# line plot
ax5 = fig1.add_subplot(gs[:, 2])
ax6 = fig1.add_subplot(gs[:, 3])


# MPI_GE
ax1, fmap, lmap = stat_overview.spatial_pattern_plot(
    ax1,
    EOFs_decade["MPI_GE"].eof.isel(decade=0),
    EOFs_decade["MPI_GE"].fra.isel(decade=0),
    EOFs_decade["MPI_GE"].eof.isel(decade=-1),
    EOFs_decade["MPI_GE"].fra.isel(decade=-1),
    levels=np.arange(-2,2.1,0.4),
)

first_eof,last_eof = split_first_last(EOFs_decade["MPI_GE"])
ax2,hist = stat_overview.index_distribution_plot(
    ax2,
    first_eof.pc,
    last_eof.pc,
)

# 20 CR 
ax3, fmap, lmap = stat_overview.spatial_pattern_plot(
    ax3,
    CR20_first_eof.eof.sel(mode="NAO").squeeze(),
    CR20_first_eof.fra.sel(mode="NAO").squeeze(),
    levels=,
)
ax4, hist = stat_overview.index_distribution_plot(
ax4,
self.first_eof_std.pc.sel(mode="NAO"),
self.last_eof_std.pc.sel(mode="NAO"),
)



ax3,lines_pos = extplt.extrc_time_line_single(
    EXTRCs,
    extr_type='pos',
    ax = ax3,
    ylim = (25,315),
    ci = True,
)

ax4,lines_neg = extplt.extrc_time_line_single(
    EXTRCs,
    extr_type='neg',
    ax = ax4,
    ylim = (25,315),
    ci = True,
)


#### ax2 ####

# add legend
f_patch = mpatches.Patch(color="#1f77b4", label="first10")
l_patch = mpatches.Patch(color="#ff7f0e", label="last10")

ax2.legend(
    handles=[f_patch, l_patch], loc="b", frameon=False,
    space = 3.2,
)

cbar = ax2.colorbar(
    fmap,
    width=0.1,
    shrink=1,
    loc = 'b', 
    label = 'standardized NAO',   
    pad = 0.2,
)

ax2.set_ylabel("density", fontsize=7)
ax2.format(
    ylocator=pplt.MultipleLocator(0.05),
    xlabel = "standard deviation",
    xlabelpad=0.8, 
    facecolor = 'none',
)


#### ax3 ####
# set the axis
ax3.spines["right"].set_visible(False)
ax3.spines["top"].set_visible(False)

ax3.format(
    xlabel="time",
    xlabelpad=0.8, 
    xtickminor=False,
    xrotation=45,
    ylim = (120,350),
    ylabel = 'Extreme occurence'
)

#### ax4 ####
# set the axis
ax4.spines["right"].set_visible(False)
ax4.spines["top"].set_visible(False)
ax4.format(
    xlabel="time",
    xlabelpad=0.8, 
    xtickminor=False,
    yticklabels = [],
    ylabel = '',
    xrotation=45,
    ylim = (120,350),
)

ax4.legend(
    lines_pos,
    models_legend,
    loc="b",
    ncols=2,
    frameon=False,
    bbox_to_anchor=(-0.1, -0.2),
    pad = 0.2,
    columnspacing=3,
)



plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/extrc_change.pdf",
    bbox_inches="tight",
)

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

axes,maps = composite_plot.plot_composite_single_ext(COMPOSITEs, models, axes)
fig3.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K",width = 0.1,shrink=1)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/composite_tsurf_pos_same_number.png",
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

axes,maps = composite_plot.plot_composite_single_ext(COMPOSITEs, models, axes,extr_type='neg')
fig4.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K",width = 0.1,shrink=1)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/composite_tsurf_neg_same_number.png",
    dpi=300,
    bbox_inches="tight",
)

# %%
