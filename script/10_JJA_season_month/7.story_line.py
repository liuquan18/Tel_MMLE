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


from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
import matplotlib.patches as mpatches


import src.composite.field_composite as field_composite
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import src.plots.utils as utils


# %%
import importlib

importlib.reload(extplt)
importlib.reload(stat_overview)
importlib.reload(field_composite)

# %%
######################
# read data
######################


def read_eof_decade(model):
    """read eofs that is decomposed by decade"""
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"
    filename = "plev_50000_decade_mpi_first_JJA_eof_result.nc"
    ds = xr.open_dataset(odir + filename)
    ds = ds.sel(mode="NAO")
    return ds


def read_eof_all(model):
    """read eofs of one model to compare with ERA5"""
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/ERA5/EOF_result/"  # ALL models are here
    filename = f"plev_50000_1940_2022_{model}_all.nc"
    ds = xr.open_dataset(odir + filename)
    ds = ds.sel(mode="NAO")
    return ds


def read_extrc(model):
    """read extreme counts"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_count/"
    filename = "plev_50000_decade_mpi_first_JJA_extre_counts.nc"
    ds = xr.open_dataset(odir + filename).pc
    return ds


def read_composite(model, var_name):
    """read composite data"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/"
    first_name = f"plev_50000_decade_mpi_first_JJA_JJA_first_{var_name}_composite.nc"
    last_name = f"plev_50000_decade_mpi_first_JJA_JJA_last_{var_name}_composite.nc"

    first = xr.open_dataset(f"{odir}{model}/composite/{first_name}")
    last = xr.open_dataset(f"{odir}{model}/composite/{last_name}")
    if var_name == "ts":
        try:
            first = first.tsurf
            last = last.tsurf
        except AttributeError:
            first = first.ts
            last = last.ts
    elif var_name == "pr":
        try:
            first = first.pr
            last = last.pr
        except AttributeError:
            first = first.precip
            last = last.precip
    return first, last


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
    elif variable == "extrc":
        all_model_data = {model: read_extrc(model) for model in models}
    elif variable == "composite":
        all_model_data = {
            model: read_composite(model, var_name="ts") for model in models
        }
    return all_model_data


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
mpl.rcParams["font.family"] = "Arial"


# %%
EOFs_decade = read_all_models("eof_decade")
EOFs_all = read_all_models("eof_all")
EXTRCs = read_all_models("extrc")
COMPOSITEs = read_all_models("composite")

# %%
# Fig 1, spatial patterns of ERA5 and MPI_GE, the time series and the boxplot
# set the fig size as 88mm x 88mm
# set the font size as 7pt
fig1 = pplt.figure(figsize=(150 / 25.4, 150 / 25.4),sharex=False,sharey=False)
fig1.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

gs = pplt.GridSpec(
    ncols=2,
    nrows=2,
    hratios=[1, 1],
    wratios = [0.6,1],
    wspace=4,
    hspace=4,
)

ax1 = fig1.add_subplot(gs[0, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax2 = fig1.add_subplot(gs[1, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax3 = fig1.add_subplot(gs[0, 1])
ax4 = fig1.add_subplot(gs[1, 1])

ax1, fmap, _ = stat_overview.spatial_pattern_plot(
    ax1,
    EOFs_all["ERA5"].eof.isel(decade=0),
    EOFs_all["ERA5"].fra.isel(decade=0),
    levels=np.arange(-40, 41, 5),
)

ax2, fmap, _ = stat_overview.spatial_pattern_plot(
    ax2,
    EOFs_all["MPI_GE"].eof.isel(decade=0),
    EOFs_all["MPI_GE"].fra.isel(decade=0),
    levels=np.arange(-40, 41, 5),
)

ax3 = stat_overview.envelop_obs_mmlea(
    ax3,
    EOFs_all["ERA5"].pc,
    EOFs_all["MPI_GE"].pc,
)

ax4 = stat_overview.obs_mmlea_box_plot(
    ax4,
    EOFs_all,
)


#### ax2 ####
cbar = ax2.colorbar(
    fmap,
    orientation="horizontal",
    shrink=1,
    ticks=np.arange(-40, 41, 10),
    ticklabelsize=6,
    labelsize=6,
    extend="both",
    loc='b',
    width = 0.1,
    pad = 1,
)
cbar.ax.set_title(
    "NAO/m",
    fontsize=7,
)




#### ax3 ####
# set the axis
ax3.spines["right"].set_visible(False)
ax3.spines["top"].set_visible(False)

# change the color line2D in ax3 into black
line_era5 = ax3.lines[0]
line_era5.set_color("black")
line_era5.set_linewidth(0.5)

# change the polycollection in ax3 into orange
poly_mpi = ax3.collections[0]
poly_mpi.set_color("orange")

# change the ticks
ax3.tick_params(
    axis="x",
    which="major",
    direction="out",
    pad=2,
    labelsize=7,
    labelcolor="black",
    labelrotation=0,
)



#### ax4 ####
# set the axis
ax4.spines["right"].set_visible(False)
ax4.spines["top"].set_visible(False)


# the ticks
ax4.tick_params(
    axis="x",
    which="major",
    direction="out",
    pad=2,
    labelsize=7,
    labelcolor="black",
    labelrotation=90,
)


plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/statistical_overview.png",
    dpi=300,
    bbox_inches="tight",
)


# %%
# Fig 2, spatial pattern chagne, index distribution change, and extreme lines
# set the fig size as 180mm x 150mm
# set the font size as 7pt
fig2 = pplt.figure(figsize=(180 / 25.4, 150 / 25.4))
fig2.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
)
models_legend = [
    "MPI_GE_onepct (100)",
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    ]

gs = pplt.GridSpec(
    ncols=3,
    nrows=2,
    wspace=(5, 0.5),
    hspace=0.5,
    hratios=[1.3, 1],
)

ax1 = fig2.add_subplot(gs[0, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax2 = fig2.add_subplot(gs[1,0])
ax3 = fig2.add_subplot(gs[:, 1])
ax4 = fig2.add_subplot(gs[:, 2])


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

# add legend
f_patch = mpatches.Patch(color="#1f77b4", label="first10")
l_patch = mpatches.Patch(color="#ff7f0e", label="last10")


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
ax2.legend(
    handles=[f_patch, l_patch], loc="b", title="periods", frameon=False,
    pad = 0.2,
)
ax2.colorbar(
    fmap,
    width=0.1,
    shrink=1,
    loc = 'b', 
    label = 'NAO_std',   
    pad = 0.2,
)

#### ax3 ####
# set the axis
ax3.spines["right"].set_visible(False)
ax3.spines["top"].set_visible(False)

#### ax4 ####
# set the axis
ax4.spines["right"].set_visible(False)
ax4.spines["top"].set_visible(False)

ax4.legend(
    lines_pos,
    models_legend,
    loc="b",
    ncols=2,
    frameon=False,
    bbox_to_anchor=(-0.1, -0.2),
    pad = 0.2,
)



plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/extrc_change.png",
    dpi=300,
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

axes,maps = field_composite.plot_composite_single_ext(COMPOSITEs, models, axes)
fig3.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K",width = 0.1,shrink=1)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/composite_tsurf_pos.png",
    dpi=300,
    bbox_inches="tight",
)

# %%
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

axes,maps = field_composite.plot_composite_single_ext(COMPOSITEs, models, axes,extr_type='neg')
fig4.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K",width = 0.1,shrink=1)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/composite_tsurf_neg.png",
    dpi=300,
    bbox_inches="tight",
)
# %%
