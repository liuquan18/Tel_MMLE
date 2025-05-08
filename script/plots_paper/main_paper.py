# %%
import xarray as xr
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import proplot as pplt
from matplotlib.lines import Line2D

from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches


import src.plots.composite_plot as composite_plot
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import src.obs.era5_extreme_change as era5_extreme_change
import src.extreme.extreme_count_troposphere as ext_profile
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

    # divide the ensemble size of each model
    ens_sizes = {
        "MPI_GE": 100,
        "MPI_GE_onepct": 100,
        "CanESM2": 50,
        "CESM1_CAM5": 40,
        "MK36": 30,
        "GFDL_CM3": 20,
    }
    ds = ds / ens_sizes[model]
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


# read slope data
def read_slope(model, x="tsurf", fixed_pattern="decade_mpi"):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_slope/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_extre_slope_{x}.nc"
    ds = xr.open_dataset(odir + filename).pc
    return ds


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
    elif variable == "slope_time":
        all_model_data = {model: read_slope(model, "time") for model in models}
    elif variable == "slope_tsurf":
        all_model_data = {model: read_slope(model, "tsurf") for model in models}
    return all_model_data


# %%
# MPI_GE random
def read_MPI_GE_random_extrc(model="MPI_GE"):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/extreme_count/"
    random_eofs = {}
    for ens_size in [20, 30, 40, 50, 80]:
        filename = odir + f"plev_50000_decade_JJA_first_{ens_size}_extre_counts.nc"
        ds = xr.open_dataset(filename)

        # divide the ensemble size
        ds = ds.pc / ens_size
        random_eofs[ens_size] = ds
    return random_eofs


def read_MPI_GE_random_slope(model="MPI_GE", x="tsurf"):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/extreme_slope/"
    random_eofs = {}
    for ens_size in [20, 30, 40, 50, 80]:
        filename = odir + f"plev_50000_decade_JJA_first_{ens_size}_extre_slope_{x}.nc"
        ds = xr.open_dataset(filename).pc

        random_eofs[ens_size] = ds
    return random_eofs


# %%
# 20CR
# eof
def read_eof_rean(model, group_size=40):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    first_eof_path = odir + "EOF_result/first_plev500_eof.nc"
    last_eof_path = odir + "EOF_result/last_plev500_eof.nc"

    first_eof = xr.open_dataset(first_eof_path)
    last_eof = xr.open_dataset(last_eof_path)
    return first_eof, last_eof


# read extreme counts
def read_extrc_rean(model, group_size=40):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    first_extc_path = odir + "first_plev50000_extc.nc"
    last_extc_path = odir + "last_plev50000_extc.nc"

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
        formater = "first40"
    elif period == 0.8:
        # formater = f"last \n 1986-2015"
        formater = "last40"
    else:
        formater = None
        print(period)
    return formater


def SMILE_rate_increase(SLOPEs, EXTRCs):
    RATEs = SLOPEs.copy()
    for key, slope in SLOPEs.items():
        try:
            first_extrc = EXTRCs[key].sel(time="1950", confidence="true")
        except KeyError:
            first_extrc = EXTRCs[key].sel(time="1950")
        rate = slope / first_extrc * 100
        RATEs[key] = rate
    return RATEs


def mean_rate(RATEs):
    # to xarray
    RATEx = []
    for key, item in RATEs.items():
        item = item.squeeze()
        RATEx.append(item)
    RATEx = xr.concat(RATEx, dim = 'model',coords='minimal',compat='override')
    # calculate the mean and std from data with mean and std
    MEANs = RATEx.sel(slopes = 'true')
    VARs = (RATEx.sel(slopes = 'high') - RATEx.sel(slopes = 'true'))**2
    MEANs = MEANs/VARs
    Weights = 1/VARs
    MEAN_best = MEANs.sum(dim = 'model')/Weights.sum(dim = 'model')

    STD_best = (Weights.sum(dim = 'model'))**0.5

    return MEAN_best, STD_best


def rean_slope(first_extrc, last_extrc, extr_type="pos", rate=True):
    first = first_extrc.sel(
        mode="NAO", confidence="true", extr_type=extr_type
    ).pc.values
    last = last_extrc.sel(mode="NAO", confidence="true", extr_type=extr_type).pc.values
    slope = (last - first) / ((1976 - 1850) / 10)

    first_high = first_extrc.sel(
        mode="NAO", confidence="high", extr_type=extr_type
    ).pc.values
    first_low = first_extrc.sel(
        mode="NAO", confidence="low", extr_type=extr_type
    ).pc.values

    last_high = last_extrc.sel(
        mode="NAO", confidence="high", extr_type=extr_type
    ).pc.values
    last_low = last_extrc.sel(
        mode="NAO", confidence="low", extr_type=extr_type
    ).pc.values

    slope_high = (last_high - first_low) / ((1976 - 1850) / 10)
    slope_low = (last_low - first_high) / ((1976 - 1850) / 10)
    if rate:
        slope = slope / first * 100
        slope_high = slope_high / first * 100
        slope_low = slope_low / first * 100

    return slope, slope_high, slope_low


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
GMST = read_all_models("gmst")

# %%
SLOPEs_time = read_all_models("slope_time", fixed_pattern=fixed_pattern)
RATEs_time = SMILE_rate_increase(SLOPEs_time, EXTRCs)
# %%
SLOPEs_tsurf = read_all_models("slope_tsurf", fixed_pattern=fixed_pattern)
RATEs_tsurf = SMILE_rate_increase(SLOPEs_tsurf, EXTRCs)

#%%
MEAN_rate_time = mean_rate(RATEs_time)
MEAN_rate_tsurf = mean_rate(RATEs_tsurf)

# %%
# MPI_GE random
MPI_GE_random_EXTRCs = read_MPI_GE_random_extrc()
MPI_GE_random_EXTRCs[100] = EXTRCs["MPI_GE"]

# MPI_GE random slope
MPI_GE_random_SLOPEs_time = read_MPI_GE_random_slope(x="time")
MPI_GE_random_SLOPEs_time[100] = SLOPEs_time["MPI_GE"]

MPI_GE_random_SLOPEs_tsurf = read_MPI_GE_random_slope(x="tsurf")
MPI_GE_random_SLOPEs_tsurf[100] = SLOPEs_tsurf["MPI_GE"]


#%%
MPI_GE_random_RATEs_time = SMILE_rate_increase(
    MPI_GE_random_SLOPEs_time, MPI_GE_random_EXTRCs
)

MPI_GE_random_RATEs_tsurf = SMILE_rate_increase(
    MPI_GE_random_SLOPEs_tsurf, MPI_GE_random_EXTRCs
)

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
CR20_first_extc, CR20_last_extc = read_extrc_rean("CR20_allens")
CR20_composite = read_composite_rean("CR20_allens", "ts")

# %%
CR20_first_extc = CR20_first_extc / (4 * 80)
CR20_last_extc = CR20_last_extc / (4 * 80)
# %% # put the 20CR_allens into the dict
COMPOSITEs["20CR"] = CR20_composite
# %%
CR20_ens_first_eof, CR20_ens_last_eof = read_eof_rean("CR20")

# also read ensemble mean of 20CR
CR20_ens_first_extc, CR20_ens_last_extc = read_extrc_rean("CR20")
CR20_ens_first_extc = CR20_ens_first_extc / 4
CR20_ens_last_extc = CR20_ens_last_extc / 4

# %%
# for 20CR
slope_pos, slope_pos_high, slope_pos_low = rean_slope(
    CR20_first_extc, CR20_last_extc, extr_type="pos", rate=False
)
slope_neg, slope_neg_high, slope_neg_low = rean_slope(
    CR20_first_extc, CR20_last_extc, extr_type="neg", rate=False
)

slope_pos_ens, slope_pos_high_ens, slope_pos_low_ens = rean_slope(
    CR20_ens_first_extc, CR20_ens_last_extc, extr_type="pos", rate=False
)
slope_neg_ens, slope_neg_high_ens, slope_neg_low_ens = rean_slope(
    CR20_ens_first_extc, CR20_ens_last_extc, extr_type="neg", rate=False
)

# %%
rate_pos, rate_pos_high, rate_pos_low = rean_slope(
    CR20_first_extc, CR20_last_extc, extr_type="pos", rate=True
)
rate_neg, rate_neg_high, rate_neg_low = rean_slope(
    CR20_first_extc, CR20_last_extc, extr_type="neg", rate=True
)

rate_pos_ens, rate_pos_high_ens, rate_pos_low_ens = rean_slope(
    CR20_ens_first_extc, CR20_ens_last_extc, extr_type="pos", rate=True
)
rate_neg_ens, rate_neg_high_ens, rate_neg_low_ens = rean_slope(
    CR20_ens_first_extc, CR20_ens_last_extc, extr_type="neg", rate=True
)

#%%
# profiles
CR20_plevs = [1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200]
CR20_mean_first_profile, CR20_mean_last_profile = ext_profile.extreme_count_rean('CR20', CR20_plevs)

onepct_first_profile, onepct_last_profile = ext_profile.extreme_count_MPI('MPI_GE_onepct')


# %%
# Fig 1, spatial pattern chagne, index distribution change, and extreme lines
# set the fig size as 180mm x 150mm
# set the font size as 7pt
fig1 = pplt.figure(figsize=(180 / 25.4, 180 / 25.4), sharex=False, sharey=False)
fig1.format(
    abc=False,
    facecolor="black",
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
    wratios=[1.3, 1.1, 0.8, 0.8],
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
    CR20_ens_first_eof.eof.sel(mode="NAO").squeeze(),
    CR20_ens_first_eof.fra.sel(mode="NAO").squeeze(),
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
    CR20_ens_first_eof.pc.sel(mode="NAO"),
    CR20_ens_last_eof.pc.sel(mode="NAO"),
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
    x=[0.3, 0.8],
    width=0.2,
    facecolor="grey",
    errcolor="grey7",
)

ax7 = extplt.reananlysis_bar(
    pos_true_first_ens,
    pos_err_first_ens,
    pos_true_last_ens,
    pos_err_last_ens,
    ax=ax7,
    x=[0.2, 0.7],
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
    x=[0.3, 0.8],
    width=0.2,
    facecolor="grey",
    errcolor="grey7",
)

ax8 = extplt.reananlysis_bar(
    neg_true_first_ens,
    neg_err_first_ens,
    neg_true_last_ens,
    neg_err_last_ens,
    ax=ax8,
    x=[0.2, 0.7],
    width=0.2,
    facecolor="none",
    errcolor="grey7",
)
ax8.set_xticks([0.2, 0.8])

ax8.xaxis.set_major_formatter(plt.FuncFormatter(format_period_year))

#### ax1 ####

cax = fig1.add_axes([0.06, 0.62, 0.23, 0.013], xlocator=MaxNLocator(5),facecolor='none')
cbar = fig1.colorbar(
    fmap_MPI,
    cax=cax,
    orientation="horizontal",
    ticks=np.arange(-30, 31, 10),
    label="Z500/m",
)
cbar.set_label("Z500/m")
cax.tick_params(axis="x", which="minor", bottom=False)
# add 'a' label at the top left corner
ax1.text(
    0.05,
    0.985,
    "a",
    transform=fig1.transFigure,
    fontweight="bold",
    va="top",
    ha="left",
)

#### ax2 ####
# add legend
ax2.axes.set_facecolor("none")
f_patch_MPI = mpatches.Patch(color="#1f77b4", label="first10")
l_patch_MPI = mpatches.Patch(color="#ff7f0e", label="last10")

ax2.set_ylabel(
    "probability density",
)
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
    fontweight="bold",
    va="top",
    ha="left",
)

first_std = first_eof.pc.std().values
last_std = last_eof.pc.std().values
ax2.set_title(
    "({first_std:0.2f} --> {last_std:0.2f})*".format(
        first_std=first_std, last_std=last_std
    ),
    pad=-15
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
    fontweight="bold",
    va="top",
    ha="left",
)


#### ax5 ####
cax = fig1.add_axes([0.06, 0.11, 0.23, 0.013], xlocator=MaxNLocator(5),facecolor='none')
cbar = fig1.colorbar(
    fmap_20CR,
    cax=cax,
    orientation="horizontal",
    ticks=np.arange(-30, 31, 10),
    label="Z500/m",
)
cbar.set_label("Z500/m")
cax.tick_params(axis="x", which="minor", bottom=False)

ax5.text(
    0.05,
    0.47,
    "e",
    transform=fig1.transFigure,
    fontweight="bold",
    va="top",
    ha="left",
)


# ax6
ax6.set_ylabel("probability density")
ax6.set_xlabel("NAO index", fontsize=7)
ax6.axes.set_facecolor("none")

ax6.legend(
    handles=[f_patch_MPI, l_patch_MPI],
    labels=["first40", "last40"],
    loc="b",
    frameon=False,
    ncol=2,
)

ax6.text(
    0.33,
    0.47,
    "f",
    transform=fig1.transFigure,
    fontweight="bold",
    va="top",
    ha="left",
)


first_std_CR20 = CR20_ens_first_eof.pc.std().values
last_std_CR20 = CR20_ens_last_eof.pc.std().values

ax6.set_title(
    "({first_std:0.2f} --> {last_std:0.2f})*".format(
        first_std=first_std_CR20, last_std=last_std_CR20
    )
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

patch_20CR = mpatches.Patch(facecolor="none", edgecolor="black", label="20CR")
patch_20CR_allens = mpatches.Patch(color="grey", label="20CR_ens (80)")

ax7.legend(
    handles=[patch_20CR, patch_20CR_allens],
    loc="b",
    frameon=False,
    ncol=2,
    bbox_to_anchor=(
        1,
        1.1,
    ),
)
ax7.text(
    0.628,
    0.47,
    "g",
    transform=fig1.transFigure,
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
    fontweight="bold",
    va="top",
    ha="left",
)


# set the axis
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)


plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_main/Fig1_MPI_GE_20CR.pdf",bbox_inches='tight'
)
# %%
# Fig 2, profiles

# %%
fig, axes = pplt.subplots(
    space=0,
    width=150 / 25.4,
    height=150 / 25.4,
    hspace=5,
    wspace=2,
    nrows=2,
    ncols=2,
    sharex=False,
    sharey=False,
)
axes.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
    xlabel="Extreme event count",
    ylabel="Pressure level (hPa)",
    ylim=(1000, 200),
    xminorticks="null",
    yminorticks="null",
    grid=False,
    toplabels=("pos", "neg"),
    # leftlabels=("NAO", "EA"),
    # xlocator=20,
)

ens_labels = ["first10", "last10"]
for i, count in enumerate([onepct_first_profile, onepct_last_profile]):
    extplt.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="pos"),
        axes[0, 0],
        label=ens_labels[i],
    )
    extplt.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="neg"),
        axes[0, 1],
        label=ens_labels[i],
    )

mean_labels = ["first40", "last40"]
for i, count in enumerate([CR20_mean_first_profile, CR20_mean_last_profile]):
    extplt.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="pos"),
        axes[1, 0],
        label=mean_labels[i],
    )
    extplt.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="neg"),
        axes[1, 1],
        label=mean_labels[i],
    )
axes[0,0].set_xlim(1,4)
axes[0,1].set_xlim(1,4)

axes[1,0].set_xlim(0,7)
axes[1,1].set_xlim(0,7)

# axes[0,0].set_xticks(np.arange(1,4.1,0.5))
# axes[0,1].set_xticks(np.arange(1,4.1,0.5))

# no yticks ylabel for the second columns
for ax in axes[:, 1]:
    ax.format(ylabel="")
plt.setp(axes[0,1].get_yticklabels(), visible=False)
plt.setp(axes[1,1].get_yticklabels(), visible=False)

axes[0, 0].legend(loc="lr", ncols=1, frame=True)
axes[1, 0].legend(loc="lr", ncols=1, frame=True)
for ax in axes:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_main/profile.pdf")


# %%
# Fig 2 line plot for other models, linear regression
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "20CR"]

fig2, axes = pplt.subplots(
    space=0,
    width=150 / 25.4,
    height=150 / 25.4,
    hspace=5,
    wspace=2,
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
    ylim=(1, 4.0),
    ci=True,
    models=["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"],
)

ax2, lines_neg_other = extplt.extrc_time_line_single(
    EXTRCs,
    extr_type="neg",
    ax=ax2,
    ylim=(1, 4.0),
    ci=True,
    models=["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"],
)

ax3, bars = extplt.extrc_slope_line(RATEs_time, ax=ax3, extr_type="pos")
ax3, bars_rand = extplt.extrc_slope_line(
    MPI_GE_random_RATEs_time,
    ax=ax3,
    extr_type="pos",
    rand=True,
)
ax4, bars = extplt.extrc_slope_line(
    RATEs_time,
    ax=ax4,
    extr_type="neg",
)

ax4, bars_rand = extplt.extrc_slope_line(
    MPI_GE_random_RATEs_time,
    ax=ax4,
    extr_type="neg",
    rand=True,
)

ax3.spines["bottom"].set_position(("data", 0))
ax4.spines["bottom"].set_position(("data", 0))


ax3.hlines(
    x1=20, x2=70, y=rate_pos, color="black", linestyle="dashed", zorder=0, linewidth=1
)
ax3.hlines(
    x1=20,
    x2=70,
    y=rate_pos_ens,
    color="black",
    linestyle="dotted",
    zorder=0,
    linewidth=1,
)

ax4.hlines(
    x1=20, x2=70, y=rate_neg, color="black", linestyle="dashed", zorder=0, linewidth=1
)
ax4.hlines(
    x1=20,
    x2=70,
    y=rate_neg_ens,
    color="black",
    linestyle="dotted",
    zorder=0,
    linewidth=1,
)


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
    ylim=(1.0, 3.8),
    xtickminor=False,
    facecolor="none",
)
ax1.set_yticks(np.arange(1.0, 3.6, 0.5))

ax2.format(
    ytickminor=False,
    grid=False,
    ylabel="",
    ylim=(1.0, 3.8),
    xtickminor=False,
    facecolor="none",
)

ax2.set_yticks(np.arange(1.0, 3.6, 0.5))

ax3.format(
    ytickminor=False,
    grid=False,
    ylabel="rate of increase per decade (%)",
    ylim=(-2.4,10.7),
    xlim = (16, 72),
    xtickminor=False,
    facecolor="none",
)
ax3.set_xticks([20, 30, 40, 50, 70])
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
    ylim=(-2.4,10.7),
    xtickminor=False,
    facecolor="none",
    xlim = (16, 72),
)
ax4.set_xticks([20, 30, 40, 50, 70])

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
    "MPI-GE (resampled)",
    "20CR_ens (80)",
    "20CR",
]
colors_model = ["red", "C1", "tab:purple", "tab:blue", "tab:green", "C4"]

legend_lines = [
    Line2D([0], [0], color=colors_model[5], lw=1.5),
    Line2D([0], [0], color=colors_model[4], lw=1.5),
    Line2D([0], [0], color=colors_model[3], lw=1.5),
    Line2D([0], [0], color=colors_model[2], lw=1.5),
    Line2D([0], [0], color=colors_model[0], lw=1.5),
    Line2D([0], [0], color=colors_model[1], lw=1.5),
    # Line2D([0], [0], color=colors_model[1], lw=1.5),
    mpatches.Patch(facecolor="none", edgecolor="C1"),
    Line2D([0], [0], color="k", lw=1.5, linestyle="dashed"),
    Line2D([0], [0], color="k", lw=1.5, linestyle="dotted"),
]


fig2.legend(
    legend_lines,
    models_legend,
    loc="b",
    frameon=False,
)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_main/Fig2_SMILEs.pdf",
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
    "20CR_ens(80)",
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
    toplabels_kw={"fontsize": 7, },
    leftlabels_kw={"fontsize": 7,},
)

axes, maps = composite_plot.plot_composite_single_ext(COMPOSITEs, models, axes)
fig3.colorbar(
    maps[0],
    loc="b",
    pad=1,
    title="(near) surface temperature / K",
    width=0.1,
    shrink=1,
)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_main/Fig3composite_pos.pdf",
    layout="tight",
)


# %%
# Fig 4, composite plot of ts for neg extremes
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "20CR"]
models_legend = [
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    "20CR_ens(80)",
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
fig4.colorbar(
    maps[0],
    loc="b",
    pad=1,
    title="(near) surface temperature / K",
    width=0.1,
    shrink=1,
)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_main/Fig4composite_neg.pdf",
    layout="tight",
)

# %%
