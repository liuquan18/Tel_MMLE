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
from matplotlib.lines import Line2D

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
#%%
# MPI_GE random
def read_MPI_GE_random_extrc(model = 'MPI_GE'):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/extreme_count/"
    random_eofs = {}
    for ens_size in [20,30,40,50]:
        filename = odir + f"plev_50000_decade_JJA_first_{ens_size}_extre_counts.nc"
        ds = xr.open_dataset(filename)
        ds = ds.pc
        random_eofs[ens_size] = ds
    return random_eofs

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
        # formater = f"first \n 1950-1979"
        formater = 'first40'
    elif period == 0.8:
        # formater = f"last \n 1986-2015"
        formater = 'last40'
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

ENS_SIZES = {'MPI_GE':100,
             'MPI_GE_onepct':100,
             'CanESM2':50,
                'CESM1_CAM5':40,
                'MK36':30,
                'GFDL_CM3':20,
                '20CR':79
             }
for model in EXTRCs.keys():
    EXTRCs[model] = EXTRCs[model]/ENS_SIZES[model]

GMST = read_all_models("gmst")

#%%
# MPI_GE random
MPI_GE_random_EXTRCs = read_MPI_GE_random_extrc()
for ens_size in MPI_GE_random_EXTRCs.keys():
    MPI_GE_random_EXTRCs[ens_size] = MPI_GE_random_EXTRCs[ens_size]/ens_size
#%%
MPI_GE_random_EXTRCs[100] = EXTRCs['MPI_GE']


# %%
# the projected spatial pattern
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
# %% # put the 20CR_allens into the dict
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
    abc=False,
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
    wspace=(4.5, 4.5, 1),
    hratios=[1, 1],
    wratios=[1.1, 1.1, 0.8, 0.8],
)

# MPI_GE
## Spatial pattern and index distribution
ax1 = fig1.add_subplot(gs[0, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax2 = fig1.add_subplot(gs[0, 1])

## line plot
ax3 = fig1.add_subplot(gs[0, 2])
ax4 = fig1.add_subplot(gs[0, 3], sharey=ax3, sharex=ax3)


# 20CR
## Spatial pattern and index distribution
ax5 = fig1.add_subplot(gs[1, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax6 = fig1.add_subplot(gs[1, 1], sharey=ax2, sharex=ax2)

## bar plot of extreme counts
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

ax5, fmap_20CR, lmap = stat_overview.spatial_pattern_plot(
    ax5,
    CR20_first_eof.eof.sel(mode="NAO").squeeze(),
    CR20_first_eof.fra.sel(mode="NAO").squeeze(),
    levels=np.arange(-30, 31, 5),
)

# index distribution
first_eof, last_eof = split_first_last(EOFs_decade["MPI_GE"])
ax2, hist = stat_overview.index_distribution_plot(
    ax2,
    first_eof.pc,
    last_eof.pc,
)
ax6, hist = stat_overview.index_distribution_plot(
    ax6,
    CR20_first_eof.pc.sel(mode="NAO"),
    CR20_last_eof.pc.sel(mode="NAO"),
)

# LINE PLOT for MPI_GE
ax3, lines_pos_MPI = extplt.extrc_time_line_single(
    EXTRCs,
    extr_type="pos",
    ax=ax3,
    ylim=(1.5, 4.5),
    ci=True,
    models=["MPI_GE_onepct", "MPI_GE"],
)

ax4, lines_neg = extplt.extrc_time_line_single(
    EXTRCs,
    extr_type="neg",
    ax=ax4,
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

cax = fig1.add_axes([0.05, 0.63, 0.23, 0.013], xlocator=MaxNLocator(5))
cbar = fig1.colorbar(
    fmap_MPI,
    cax=cax,
    orientation="horizontal",
    ticks=np.arange(-30, 31, 10),
    label="Z500/m",
)
cax.tick_params(axis='x', which='minor', bottom=False)
# add 'a' label at the top left corner
ax1.text(
    0.05,
    0.985,
    "a",
    transform=fig1.transFigure,
    fontsize=7,
    fontweight="bold",
    va="top",
    ha="left",
)

#### ax2 ####
# add legend
ax2.axes.set_facecolor("none")
f_patch_MPI = mpatches.Patch(color="#1f77b4", label="first10")
l_patch_MPI = mpatches.Patch(color="#ff7f0e", label="last10")

ax2.set_ylabel("probability density", fontsize=7,)
ax2.set_xlabel("NAO index", fontsize=7)

ax2.legend(
    handles=[f_patch_MPI, l_patch_MPI],
    loc="b",
    frameon=False,
    ncol=2,
)
# add 'b' to ax2
ax2.text(
    0.33,
    0.985,
    "b",
    transform=fig1.transFigure,
    fontsize=7,
    fontweight="bold",
    va="top",
    ha="left",
)

#### ax3 ####
# set the axis
ax3.spines["right"].set_visible(False)
ax3.spines["top"].set_visible(False)

ax3.format(
    xlabel="",
    xlabelpad=0.8,
    xtickminor=False,
    xrotation=45,
    ylabel="Extreme occurence decade$^{-1}$ realization$^{-1}$",
    ylim=(1.5, 3.6),
    facecolor="none",
)
ax3.set_yticks(np.arange(1.5, 3.5, 0.5))


ax3.legend(
    lines_pos_MPI,
    models_legend,
    loc="b",
    ncols=1,
    bbox_to_anchor=(
        1.2,
        1.15,
    ),
    frameon=False,
    # major tick every 0.5
)

ax3.text(
    0.628,
    0.985,
    "c",
    transform=fig1.transFigure,
    fontsize=7,
    fontweight="bold",
    va="top",
    ha="left",
)

#### ax4 ####
# set the axis
ax4.spines["right"].set_visible(False)
ax4.spines["top"].set_visible(False)
ax4.format(
    xlabel="",
    xlabelpad=0.8,
    xtickminor=False,
    ylabel="",
    ylim=(1.5, 3.6),
    xrotation=45,
    facecolor="none",
)

ax4.text(
    0.81,
    0.985,
    "d",
    transform=fig1.transFigure,
    fontsize=7,
    fontweight="bold",
    va="top",
    ha="left",
)

#### ax5 ####
cax = fig1.add_axes([0.05, 0.13, 0.23, 0.013], xlocator=MaxNLocator(5))
cbar = fig1.colorbar(
    fmap_20CR,
    cax=cax,
    orientation="horizontal",
    ticks=np.arange(-30, 31, 10),
    label="Z500/m",
)
cax.tick_params(axis='x', which='minor', bottom=False)

ax5.text(
    0.05,
    0.47,
    "e",
    transform=fig1.transFigure,
    fontsize=7,
    fontweight="bold",
    va="top",
    ha="left",
)


# ax6
ax6.set_ylabel("probability density", fontsize=7)
ax6.set_xlabel("NAO index", fontsize=7)
ax6.axes.set_facecolor("none")

ax6.legend(
    handles=[f_patch_MPI, l_patch_MPI],
    labels = ["first40", "last40"],
    loc="b",
    frameon=False,
    ncol=2,
)

ax6.text(
    0.33,
    0.47,
    "f",
    transform=fig1.transFigure,
    fontsize=7,
    fontweight="bold",
    va="top",
    ha="left",
)

#### ax7 ####""
ax7.format(
    ytickminor=False,
    grid=False,
    ylabel="Extreme occurence decade$^{-1}$ realization$^{-1}$",
    ylim=(0, 6.5),
    xlim=(0, 1),
    xtickminor=False,
    facecolor="none",
)
ax7.spines["right"].set_visible(False)
ax7.spines["top"].set_visible(False)

patch_20CR = mpatches.Patch(facecolor="none",edgecolor='black', label="20CR")
patch_20CR_allens = mpatches.Patch(color="grey", label="20CR_allens (79)")

ax7.legend(
    handles=[patch_20CR,patch_20CR_allens],
    loc="b",
    frameon=False,
    ncol=2,
    bbox_to_anchor=(
        1.2,
        1.1,
    ),
)
ax7.text(
    0.628,
    0.47,
    "g",
    transform=fig1.transFigure,
    fontsize=7,
    fontweight="bold",
    va="top",
    ha="left",
)

#### ax8 ####
ax8.format(
    ytickminor=False,
    grid=False,
    ylabel="",
    ylim=(0, 6.5),
    xtickminor=False,
    xlim=(0, 1),
    facecolor="none",
)
ax8.spines["right"].set_visible(False)
ax8.spines["top"].set_visible(False)

ax8.text(
    0.81,
    0.47,
    "h",
    transform=fig1.transFigure,
    fontsize=7,
    fontweight="bold",
    va="top",
    ha="left",
)

# set the axis
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)


# plt.savefig(
#     "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/Fig1_MPI_GE_20CR.pdf",
# )
# %%
# Fig 2 line plot for other models, linear regression
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "20CR"]

fig2, axes = pplt.subplots(
    space=0,
    width=150 / 25.4,
    height=150 / 25.4,
    hspace = 5,
    wspace = 2,
    nrows=2,
    ncols=2,
    sharex=False,
    sharey=False,
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
    ylim=(1,4.0),
    ci=True,
    models=["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"],
)

ax2, lines_neg_other = extplt.extrc_time_line_single(
    EXTRCs,
    extr_type="neg",
    ax=ax2,
    ylim=(1,4.0),
    ci=True,
    models=["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"],
)

ax3,bars = extplt.extrc_slope_line(EXTRCs,ax = ax3,tsurfs = None,extr_type ='pos',time = '1960-06-06')
ax3,bars_rand = extplt.extrc_slope_line(MPI_GE_random_EXTRCs,ax = ax3,tsurfs = None,extr_type ='pos',rand = True,time = '1960-06-06')
ax4,bars = extplt.extrc_slope_line(EXTRCs,ax = ax4,tsurfs = None,extr_type ='neg',time = '1960-06-06')
ax4,bars_rand = extplt.extrc_slope_line(MPI_GE_random_EXTRCs,ax = ax4,tsurfs = None,extr_type ='neg',rand = True,time = '1960-06-06')

# hline at y = 0 for ax3 and ax4
# ax3.axhline(y=0, color="black", linewidth = ax3.spines['bottom'].get_linewidth())
# ax4.axhline(y=0, color="black", linewidth = ax4.spines['bottom'].get_linewidth())

ax3.spines['bottom'].set_position(('data', 0))
ax4.spines['bottom'].set_position(('data', 0))

###
for ax in [ax1, ax2, ax3, ax4]:
    # no right and top axis
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

### ax1
ax1.format(
    ytickminor=False,
    grid=False,
    ylabel="Extreme occurence decade$^{-1}$ realization$^{-1}$",
    ylim=(1.0,3.8),
    xtickminor=False,
    facecolor="none",
)
ax1.set_yticks(np.arange(1.0, 3.6, 0.5))

ax2.format(
    ytickminor=False,
    grid=False,
    ylabel="",
    ylim=(1.0,3.8),
    xtickminor=False,
    facecolor="none",
)

ax2.set_yticks(np.arange(1.0, 3.6, 0.5))

ax3.format(
    ytickminor=False,
    grid=False,
    ylabel="increase in occurence decade$^{-1}$ realization$^{-1}$",
    ylim=(-0.049,0.18),
    xtickminor=False,
    facecolor="none",
)
ax3.set_xticks([20,30,40,50,70])
# the xlabel for ax3 at the right bottom
ax3.text(
    0.98,
    0.02,
    "ensemble size",
    transform=ax3.transAxes,
    va="bottom",
    ha="right",
)


ax4.format(
    ytickminor=False,
    grid=False,
    ylabel="",
    ylim=(-0.049,0.18),
    xtickminor=False,
    facecolor="none",
)
ax4.set_xticks([20,30,40,50,70])
ax4.text(
    0.98,
    0.02,
    "ensemble size",
    transform=ax4.transAxes,
    va="bottom",
    ha="right",
)


plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)

models_legend = [
    "GFDL-CM3 (20)",
    "MK3.6 (30)",
    "CESM1-CAM5 (40)",
    "CanESM2 (50)",
    "MPI-GE_onepct (100)",
    "MPI-GE (100)",
    "MPI-GE resampled"
]
colors_model = ["red", "C1", "tab:purple", "tab:blue", "tab:green", "C4"]

legend_lines = [
    Line2D([0], [0], color=colors_model[5], lw=1.5),
    Line2D([0], [0], color=colors_model[4], lw=1.5),
    Line2D([0], [0], color=colors_model[3], lw=1.5),
    Line2D([0], [0], color=colors_model[2], lw=1.5),
    Line2D([0], [0], color=colors_model[0], lw=1.5),
    Line2D([0], [0], color=colors_model[1], lw=1.5),
    Line2D([0], [0], color=colors_model[1], lw=1.5,linestyle='--'),
]


fig2.legend(
    legend_lines,
    models_legend,
    loc="b",
    frameon=False,
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
