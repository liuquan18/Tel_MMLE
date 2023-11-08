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
def read_eof_decade(model, fixed_pattern="decade_mpi"):
    """read eofs that is decomposed by decade"""
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_eof_result.nc"
    ds = xr.open_dataset(odir + filename)
    ds = ds.sel(mode="NAO")
    return ds


def read_extrc(model, fixed_pattern="decade_mpi"):
    """read extreme counts"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_count/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_extre_counts.nc"
    ds = xr.open_dataset(odir + filename).pc
    return ds


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


def read_gmst(model, tsurf="ens_fld_year_mean"):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    ts = xr.open_dataset(f"{odir}{tsurf}.nc")
    try:
        ts = ts.tsurf
    except AttributeError:
        try:
            ts = ts.ts
        except AttributeError:
            ts = ts.tas
    return ts


def read_all_models(variable, fixed_pattern="decade_mpi"):
    """read all models data"""
    models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    if variable == "eof_decade":
        all_model_data = {
            model: read_eof_decade(model, fixed_pattern) for model in models
        }
    elif variable == "extrc":
        all_model_data = {model: read_extrc(model, fixed_pattern) for model in models}
    elif variable == "composite":
        models = models[1:]  # no MPI_GE_onepct
        all_model_data = {
            model: read_composite(model, var_name="ts") for model in models
        }
    elif variable == "gmst":
        all_model_data = {model: read_gmst(model) for model in models}
    return all_model_data


# %%
# 20CR
# eof
def read_eof_rean(model, group_size=40):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    first_eof_path = odir + f"EOF_result/first_{str(group_size)}_eof_std.nc"
    last_eof_path = odir + f"EOF_result/last_{str(group_size)}_eof_std.nc"

    first_eof = xr.open_dataset(first_eof_path)
    last_eof = xr.open_dataset(last_eof_path)
    return first_eof, last_eof


# read extreme counts
def read_extrc_rean(model, group_size=40):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    first_extc_path = odir + f"first_{str(group_size)}_extc.nc"
    last_extc_path = odir + f"last_{str(group_size)}_extc.nc"

    first_extc = xr.open_dataset(first_extc_path)
    last_extc = xr.open_dataset(last_extc_path)
    return first_extc, last_extc


# read composite
def read_composite_rean(model, var_name, reduction="mean", group_size=40):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/composite/"
    composite_path = odir + "composite_mean_ts_40_withboot.nc"
    composite = xr.open_dataset(composite_path)
    return composite.__xarray_dataarray_variable__


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


def y_yerr(extrc):
    pos_true = extrc.sel(mode="NAO", confidence="true", extr_type="pos").pc.values
    pos_high = extrc.sel(mode="NAO", confidence="high", extr_type="pos").pc.values
    pos_low = extrc.sel(mode="NAO", confidence="low", extr_type="pos").pc.values
    pos_err = [pos_true - pos_low, pos_high - pos_true]

    neg_true = extrc.sel(mode="NAO", confidence="true", extr_type="neg").pc.values
    neg_high = extrc.sel(mode="NAO", confidence="high", extr_type="neg").pc.values
    neg_low = extrc.sel(mode="NAO", confidence="low", extr_type="neg").pc.values
    neg_err = [neg_true - neg_low, neg_high - neg_true]

    return pos_true, pos_err, neg_true, neg_err


def format_period_year(period, tick_number):
    if period == 0.2:
        formater = f"first \n 1950-1979"
    elif period == 0.8:
        formater = f"last \n 1986-2015"
    else:
        formater = None
        print(period)
    return formater


# %%
######################
# matplotlib global parameters
######################

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
plt.rcParams["hatch.linewidth"] = 0.3
# %%
# SMILEs
fixed_pattern = "decade_mpi"
EOFs_decade = read_all_models("eof_decade", fixed_pattern=fixed_pattern)
EXTRCs = read_all_models("extrc", fixed_pattern=fixed_pattern)
COMPOSITEs = read_all_models("composite", fixed_pattern=fixed_pattern)

# divided by 100 for each model in EXTRCs
for model in EXTRCs.keys():
    EXTRCs[model] = EXTRCs[model] / 100

# %%
GMST = read_all_models("gmst")
# %%
first_pattern = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/first_pattern_projected.nc"
).__xarray_dataarray_variable__
last_pattern = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/last_pattern_projected.nc"
).__xarray_dataarray_variable__

# %%
# 20CR all ens
CR20_first_eof, CR20_last_eof = read_eof_rean("CR20_allens")
CR20_first_extc, CR20_last_extc = read_extrc_rean("CR20_allens")
CR20_composite = read_composite_rean("CR20_allens", "ts")

# %%
CR20_first_extc = CR20_first_extc / (4 * 79)
CR20_last_extc = CR20_last_extc / (4 * 79)
# %%
COMPOSITEs["20CR"] = CR20_composite
# %%
# also read ensemble mean of 20CR
CR20_ens_first_extc, CR20_ens_last_extc = read_extrc_rean("CR20")
# %%
CR20_ens_first_extc = CR20_ens_first_extc / 4
CR20_ens_last_extc = CR20_ens_last_extc / 4
# %%
# Fig 1, spatial pattern chagne, index distribution change, and extreme lines
# set the fig size as 180mm x 150mm
# set the font size as 7pt
fig1 = pplt.figure(figsize=(180 / 25.4, 180 / 25.4), sharex=False, sharey=False)
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
    hspace=3.5,
    wspace=(3, 4.5, 1),
    hratios=[1, 1],
    wratios=[1.1, 1.1, 0.8, 0.8],
)

# spatial pattern
ax1 = fig1.add_subplot(gs[0, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax3 = fig1.add_subplot(gs[0, 1])

# index distribution
ax2 = fig1.add_subplot(gs[1, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax4 = fig1.add_subplot(gs[1, 1], sharey=ax3, sharex=ax3)

# line plot for MPI_GE
ax5 = fig1.add_subplot(gs[0, 2])
ax6 = fig1.add_subplot(gs[0, 3], sharey=ax5, sharex=ax5)

# line plot for others
ax7 = fig1.add_subplot(gs[1, 2])
ax8 = fig1.add_subplot(gs[1, 3], sharey=ax7, sharex=ax7)

# spatial pattern
ax1, fmap_MPI, lmap = stat_overview.spatial_pattern_plot(
    ax1,
    first_pattern,
    EOFs_decade["MPI_GE"].fra.isel(decade=0),
    last_pattern,
    EOFs_decade["MPI_GE"].fra.isel(decade=-1),
    levels=np.arange(-30, 31, 5),
)

ax2, fmap_20CR, lmap = stat_overview.spatial_pattern_plot(
    ax2,
    CR20_first_eof.eof.sel(mode="NAO").squeeze(),
    CR20_first_eof.fra.sel(mode="NAO").squeeze(),
    levels=np.arange(-30, 31, 5),
)

# index distribution
first_eof, last_eof = split_first_last(EOFs_decade["MPI_GE"])
ax3, hist = stat_overview.index_distribution_plot(
    ax3,
    first_eof.pc,
    last_eof.pc,
)
ax4, hist = stat_overview.index_distribution_plot(
    ax4,
    CR20_first_eof.pc.sel(mode="NAO"),
    CR20_last_eof.pc.sel(mode="NAO"),
)

# LINE PLOT for MPI_GE
ax5, lines_pos_MPI = extplt.extrc_time_line_single(
    EXTRCs,
    extr_type="pos",
    ax=ax5,
    ylim=(1.5, 4.5),
    ci=True,
    models=["MPI_GE_onepct", "MPI_GE"],
)

ax6, lines_neg = extplt.extrc_time_line_single(
    EXTRCs,
    extr_type="neg",
    ax=ax6,
    ylim=(1.5, 4.5),
    ci=True,
    models=["MPI_GE_onepct", "MPI_GE"],
)

# error bar for 20CR
pos_true_first, pos_err_first, neg_true_first, neg_err_first = y_yerr(CR20_first_extc)
pos_true_last, pos_err_last, neg_true_last, neg_err_last = y_yerr(CR20_last_extc)

pos_true_first_ens, pos_err_first_ens, neg_true_first_ens, neg_err_first_ens = y_yerr(
    CR20_ens_first_extc
)
pos_true_last_ens, pos_err_last_ens, neg_true_last_ens, neg_err_last_ens = y_yerr(
    CR20_ens_last_extc
)

ax7 = extplt.reananlysis_bar(
    pos_true_first,
    pos_err_first,
    pos_true_last,
    pos_err_last,
    ax=ax7,
    x=[0.2, 0.7],
    width=0.2,
    facecolor="grey",
)

ax7 = extplt.reananlysis_bar(
    pos_true_first_ens,
    pos_err_first_ens,
    pos_true_last_ens,
    pos_err_last_ens,
    ax=ax7,
    x=[0.3, 0.8],
    width=0.2,
    facecolor="none",
    errcolor="grey7",
)

ax7.xaxis.set_major_formatter(plt.FuncFormatter(format_period_year))
ax7.set_xticks([0.2, 0.8])

ax8 = extplt.reananlysis_bar(
    neg_true_first,
    neg_err_first,
    neg_true_last,
    neg_err_last,
    ax=ax8,
    x=[0.2, 0.7],
    width=0.2,
    facecolor="grey",
)

ax8 = extplt.reananlysis_bar(
    neg_true_first_ens,
    neg_err_first_ens,
    neg_true_last_ens,
    neg_err_last_ens,
    ax=ax8,
    x=[0.3, 0.8],
    width=0.2,
    facecolor="none",
    errcolor="grey7",
)
ax8.set_xticks([0.2, 0.8])


ax8.xaxis.set_major_formatter(plt.FuncFormatter(format_period_year))


#### ax1 ####
cbar = ax1.colorbar(fmap_MPI, width=0.1, shrink=1, loc="b", pad=0)
cbar.set_label("Z500/m", labelpad=-40, y=1.05, rotation=0)

#### ax2 ####
cbar = ax2.colorbar(
    fmap_20CR,
    width=0.1,
    shrink=1,
    loc="b",
    pad=0,
)

cbar.set_label("Z500/m", labelpad=-40, y=1.05, rotation=0)

#### ax3 ####
# add legend
f_patch = mpatches.Patch(color="#1f77b4", label="first")
l_patch = mpatches.Patch(color="#ff7f0e", label="last")


ax3.set_ylabel("density", fontsize=8, loc="top", labelpad=-15, rotation=0)
ax3.set_xlabel("NAO index", fontsize=7)

# ax4
ax4.set_ylabel("density", fontsize=8, loc="top", labelpad=-15, rotation=0)
ax4.set_xlabel("NAO index", fontsize=7)

ax4.legend(
    handles=[f_patch, l_patch],
    loc="b",
    frameon=False,
    bbox_to_anchor=(
        0.3,
        -0.3,
    ),
    ncol=2,
)

#### ax5 ####
# set the axis
ax5.spines["right"].set_visible(False)
ax5.spines["top"].set_visible(False)

ax5.format(
    xlabel="",
    xlabelpad=0.8,
    xtickminor=False,
    xrotation=45,
    ylabel="Extreme occurence / 10 years",
    ylim=(1.5, 3.5),
)
ax5.yaxis.set_major_locator(ticker.AutoLocator())

#### ax4 ####
# set the axis
ax6.spines["right"].set_visible(False)
ax6.spines["top"].set_visible(False)
ax6.format(
    xlabel="",
    xlabelpad=0.8,
    xtickminor=False,
    ylabel="",
    ylim=(1.5, 3.5),
    xrotation=45,
)
ax6.yaxis.set_major_locator(ticker.AutoLocator())


#### ax7 ####
ax7.format(
    ytickminor=False,
    grid=False,
    ylabel="Extreme occurence / 10 years",
    ylim=(0, 6.5),
    xlim=(0, 1),
    xtickminor=False,
)
ax7.spines["right"].set_visible(False)
ax7.spines["top"].set_visible(False)
ax7.legend(
    lines_pos_MPI[0],
    models_legend[0],
    loc="b",
    ncols=1,
    frameon=False,
    bbox_to_anchor=(
        -0.2,
        -0.3,
    ),
)

#### ax8 ####
ax8.format(
    ytickminor=False,
    grid=False,
    ylabel="",
    ylim=(0, 6.5),
    xtickminor=False,
    xlim=(0, 1),
)
ax8.spines["right"].set_visible(False)
ax8.spines["top"].set_visible(False)
ax8.legend(
    lines_pos_MPI[1],
    models_legend[1],
    loc="b",
    ncols=1,
    frameon=False,
    bbox_to_anchor=(
        0.3,
        -0.3,
    ),
)
# set the axis
plt.setp(ax6.get_yticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)


# plt.savefig(
#     "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/extrc_change.pdf",
#     bbox_inches="tight",
# )
# %%
# Fig 2 line plot for other models, linear regression
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "20CR"]
models_legend = [
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    "20CR(79)",
]

fig2, axes = pplt.subplots(
    space=0,
    width=180 / 25.4,
    height=180 / 25.4,
    wspace=0.2,
    hspace=0.2,
    nrows=2,
    ncols=3,
)

fig2.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

ax1 = axes[0, 0]
ax2 = axes[0, 1]

ax3 = axes[1, 0]
ax4 = axes[1, 1]


ax1, lines_pos_other = extplt.extrc_time_line_single(
    EXTRCs,
    extr_type="pos",
    ax=ax1,
    ylim=(0, 2),
    ci=True,
    models=["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"],
)

ax2, lines_neg_other = extplt.extrc_time_line_single(
    EXTRCs,
    extr_type="neg",
    ax=ax2,
    ylim=(0, 2.5),
    ci=True,
    models=["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"],
)

ext_plot.mmle_slope_scatter(
    extrs, tsurfs, extrs_rand, tsurfs_rand, tsurf=tsurf, time=time
)


###
ax1.format(
    ytickminor=False,
    grid=False,
    ylabel="Extreme occurence / 10 years",
    ylim=(0, 1.5),
    xtickminor=False,
)

ax2.format(
    ytickminor=False,
    grid=False,
    ylabel="",
    ylim=(0, 1.5),
    xtickminor=False,
)


# %%
# Fig 3, composite plot of ts for positve extremes
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "20CR"]
models_legend = [
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    "20CR(79)",
]

fig3, axes = pplt.subplots(
    space=0,
    width=180 / 25.4,
    wspace=0.2,
    hspace=0.2,
    proj="ortho",
    proj_kw=({"lon_0": -20, "lat_0": 60}),
    nrows=3,
    ncols=6,
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
    toplabels_kw={"fontsize": 7},
    leftlabels_kw={"fontsize": 7},
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

axes, maps = composite_plot.plot_composite_single_ext(COMPOSITEs, models, axes)
fig3.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K", width=0.1, shrink=1)

# plt.savefig(
#     "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/composite_tsurf_pos_same_number.png",
#     dpi=300,
#     bbox_inches="tight",
# )


# %%
# Fig 4, composite plot of ts for neg extremes
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "20CR"]
models_legend = [
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    "20CR(79)",
]

fig4, axes = pplt.subplots(
    space=0,
    width=180 / 25.4,
    wspace=0.2,
    hspace=0.2,
    proj="ortho",
    proj_kw=({"lon_0": -20, "lat_0": 60}),
    nrows=3,
    ncols=6,
)
axes.format(
    latlines=20,
    lonlines=30,
    color="grey7",
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
    leftlabels=["first10", "last10", "last10 - first10"],
    toplabels=models_legend,
    toplabels_kw={"fontsize": 7},
    leftlabels_kw={"fontsize": 7},
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

axes, maps = composite_plot.plot_composite_single_ext(
    COMPOSITEs, models, axes, extr_type="neg"
)
fig4.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K", width=0.1, shrink=1)

# plt.savefig(
#     "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/composite_tsurf_neg_same_number.png",
#     dpi=300,
#     bbox_inches="tight",
# )

# %%
