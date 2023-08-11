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
mpl.rcParams['font.family'] = 'sans-serif'

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
    wspace=5,
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

ax3 = era5_extreme_change.plot_era_nao_index(
    EOFs_all["ERA5"].pc,
    EOFs_all["ERA5_no_dec"].pc,
    ax3,
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
    label = '',
)
cbar.ax.set_title(
    "NAO/m",
    fontsize=7,
)




#### ax3 ####
# set the axis
ax3.spines["right"].set_visible(False)
ax3.spines["top"].set_visible(False)

# change the ticks
ax3.tick_params(
    axis="x",
    which="major",
    direction="out",
    pad=2,
    labelsize=7,
    labelcolor="black",
)
ax3.format(
    ylabel = 'std_NAO',
    grid = False,
    yminorticks="null",
    xminorticks="null",
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
    labelrotation=45,
)
ax4.format(
    ylabel = 'std_NAO',
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
fig2 = pplt.figure(figsize=(180 / 25.4, 180 / 25.4),sharex=False,sharey=False)
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
    hratios=[1, 0.8],
    wratios=[1, 0.85, 0.85],
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
    label = 'NAO_std',   
    pad = 0.2,
)

ax2.set_ylabel("density", fontsize=7)
ax2.format(
    ylocator=pplt.MultipleLocator(0.05),
    xlabel = "std_NAO",
    xlabelpad=0.8, 
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

axes,maps = composite_plot.plot_composite_single_ext(COMPOSITEs, models, axes)
fig3.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K",width = 0.1,shrink=1)

# plt.savefig(
#     "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/composite_tsurf_pos_same_number.png",
#     dpi=300,
#     bbox_inches="tight",
# )

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

axes,maps = composite_plot.plot_composite_single_ext(COMPOSITEs, models, axes,extr_type='neg')
fig4.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K",width = 0.1,shrink=1)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/composite_tsurf_neg_same_number.png",
    dpi=300,
    bbox_inches="tight",
)
# %%
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# create a sample data array with values of 1 and 0
data = np.random.randint(0, 2, size=(10, 10))
da = xr.DataArray(data)

# create a contour plot with hatching on the significant areas
fig, axes = pplt.subplots(
    space=0,
    width = 180/25.4,
    wspace=0.2,
    hspace=0.2,
    proj="ortho",
    proj_kw=({"lon_0": -20, "lat_0": 60}),
    nrows=3,
    ncols=1,
)
model = 'MPI_GE'
extr_type = 'pos'
first = COMPOSITEs[model].sel(mode='NAO',period = 'first', extr_type = extr_type)
last = COMPOSITEs[model].sel(mode='NAO',period = 'last', extr_type = extr_type)
diff = COMPOSITEs[model].sel(mode='NAO',period = 'diff', extr_type = extr_type)
diff_sig = COMPOSITEs[model].sel(mode='NAO',period = 'diff_sig', extr_type = extr_type)

first = utils.erase_white_line(first)
last = utils.erase_white_line(last)
diff = utils.erase_white_line(diff)
diff_sig = utils.erase_white_line(diff_sig)

data_all = [
first,
last,
diff,
]
maps = []
for j, data in enumerate(data_all):
    cs = axes[j].contourf(data, levels = np.arange(-1.5,1.6,0.3),cmap="RdBu_r",)
    axes[j].set_aspect('equal')

axes[2].contourf(diff_sig, levels=[-0.5, 0.5, 1.5], colors=['none', 'none'], hatches=['', '///'])

plt.show()
# %%
