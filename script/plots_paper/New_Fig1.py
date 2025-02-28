#%%
import xarray as xr
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import proplot as pplt
from matplotlib.lines import Line2D
import cartopy.crs as ccrs

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
# %%

def read_all_models(variable, fixed_pattern="decade_mpi"):
    """read all models data"""
    models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    if variable == "eof_decade":
        all_model_data = {
            model: read_eof_decade(model, fixed_pattern) for model in models
        }
    elif variable == "extrc":
        all_model_data = {model: read_extrc(model, fixed_pattern) for model in models}
    return all_model_data

# %%
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
# %%
# SMILEs
fixed_pattern = "decade_mpi"
EOFs_decade = read_all_models("eof_decade", fixed_pattern=fixed_pattern)
EXTRCs = read_all_models("extrc", fixed_pattern=fixed_pattern)
# %%
# 20CR all ens
CR20_first_extc, CR20_last_extc = read_extrc_rean("CR20_allens")

CR20_first_extc = CR20_first_extc / (4 * 80)
CR20_last_extc = CR20_last_extc / (4 * 80)

# %%
CR20_ens_first_eof, CR20_ens_last_eof = read_eof_rean("CR20")

CR20_ens_first_pc = CR20_ens_first_eof.pc.sel(mode="NAO") - CR20_ens_first_eof.pc.sel(mode = 'NAO').mean()
CR20_ens_last_pc = CR20_ens_last_eof.pc.sel(mode="NAO") - CR20_ens_last_eof.pc.sel(mode = 'NAO').mean()

# also read ensemble mean of 20CR
CR20_ens_first_extc, CR20_ens_last_extc = read_extrc_rean("CR20")
CR20_ens_first_extc = CR20_ens_first_extc / 4
CR20_ens_last_extc = CR20_ens_last_extc / 4

# %%
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
    ncols=3,
    nrows=2,
    hspace=3.5,
    wspace=( 4.5, 1),
    hratios=[1, 1],
    wratios=[1.1, 0.8, 0.8],
)

# MPI_GE
## Spatial pattern and index distribution
ax2 = fig1.add_subplot(gs[0, 0])

## line plot
ax3 = fig1.add_subplot(gs[0, 1])
ax4 = fig1.add_subplot(gs[0, 2], sharey=ax3, sharex=ax3)


# 20CR
## Spatial pattern and index distribution
ax6 = fig1.add_subplot(gs[1, 0], sharey=ax2, sharex=ax2)

## bar plot of extreme counts
ax7 = fig1.add_subplot(gs[1, 1])
ax8 = fig1.add_subplot(gs[1, 2], sharey=ax7, sharex=ax7)

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



# set the axis
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)

# add a, b, c
ax2.text(-0.1, 1.05, "a", transform=ax2.transAxes, fontsize=12, fontweight="bold")
ax3.text(-0.1, 1.05, "b", transform=ax3.transAxes, fontsize=12, fontweight="bold")
ax4.text(-0.1, 1.05, "c", transform=ax4.transAxes, fontsize=12, fontweight="bold")

ax6.text(-0.1, 1.05, "d", transform=ax6.transAxes, fontsize=12, fontweight="bold")
ax7.text(-0.1, 1.05, "e", transform=ax7.transAxes, fontsize=12, fontweight="bold")
ax8.text(-0.1, 1.05, "f", transform=ax8.transAxes, fontsize=12, fontweight="bold")

plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_main/New_Fig1_index.pdf", dpi = 300)

# %%
