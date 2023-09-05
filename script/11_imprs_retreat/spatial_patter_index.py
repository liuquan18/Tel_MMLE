# %%
import pandas as pd
import matplotlib as mpl
import proplot as pplt
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import matplotlib.ticker as ticker

# %%
import src.plots.statistical_overview as stat_overview
import src.plots.extreme_plot as extplt
import src.plots.composite_plot as composite_plot


#%%
import importlib
importlib.reload(stat_overview)
importlib.reload(extplt)
#%%


mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
plt.rcParams["hatch.linewidth"] = 0.3

# black background
mpl.rcParams["axes.facecolor"] = "black"
mpl.rcParams["figure.facecolor"] = "black"
mpl.rcParams["savefig.facecolor"] = "black"

# white font color
mpl.rcParams["text.color"] = "white"
mpl.rcParams["axes.labelcolor"] = "white"
mpl.rcParams["axes.edgecolor"] = "white"
mpl.rcParams["xtick.color"] = "white"
mpl.rcParams["ytick.color"] = "white"
plt.rcParams["xtick.color"] = "white"
plt.rcParams["ytick.color"] = "white"

# %%
############################
# spatial pattern and index distribution
EOFs = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result.nc"
)
EOFs = EOFs.sel(mode="NAO")
#%%

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
fig2 = pplt.figure(figsize=(180 / 25.4, 90 / 25.4),sharex=False,sharey=False)

params = {
    "figure.facecolor": "black",
    "ytick.color": "w",
    "xtick.color": "w",
    "axes.labelcolor": "w",
    "axes.edgecolor": "w",
    "tick.labelcolor": "w",
    "text.color": "w",}
pplt.rc.update(params)


gs = pplt.GridSpec(
    ncols=2,
    nrows=1,
)

ax1 = fig2.add_subplot(gs[0, 0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
ax2 = fig2.add_subplot(gs[0,1])


ax1, fmap, lmap = stat_overview.spatial_pattern_plot(
    ax1,
    EOFs.eof.isel(decade=0),
    EOFs.fra.isel(decade=0),
    # EOFs.eof.isel(decade=-1),
    # EOFs.fra.isel(decade=-1),
    levels=np.arange(-2,2.1,0.4),
)
ax1.set_facecolor("white")


first_eof,last_eof = split_first_last(EOFs)
ax2,hist = stat_overview.index_distribution_plot(
    ax2,
    first_eof.pc,
    # last_eof.pc,
)

ax2.set_facecolor("black")
ax2.set_xlabel("")
ax2.set_ylabel("probability")

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/MPI_GE_stats_first.pdf",
    bbox_inches="tight",
    dpi=300,
)


# %%
# %%
############################
# extreme count evolution
extrc =xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/extreme_count/plev_50000_decade_mpi_first_JJA_extre_counts.nc"
)
extrc = extrc.pc


#%%
fig = pplt.figure(figsize=(180 / 25.4, 90 / 25.4),sharex=False,sharey=False)

fig.format(
    tight=False,
)


params = {
    "figure.facecolor": "black",
    # "axes.facecolor": "black",
    "ytick.color": "w",
    "xtick.color": "w",
    "axes.labelcolor": "w",
    "axes.edgecolor": "w",
    "tick.labelcolor": "w",
    "text.color": "w",}
pplt.rc.update(params)


gs = pplt.GridSpec(
    ncols=2,
    nrows=1,
)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0,1])

def plot_extreme_line(extrc,ax,extr_type='pos',mode = 'NAO',ci = True):
    line = extrc.sel(extr_type=extr_type,mode = mode,confidence = 'true').plot.line(
        ax=ax, 
        label='MPI-GE',
        x = 'time',
        color = 'C1',
        linewidth = 2,)

    if ci:
    # fill between the confidence interval ['low','high']
        ax.fill_between(
        extrc.time,
        extrc.sel(extr_type=extr_type,mode = mode,confidence = 'low').values,
        extrc.sel(extr_type=extr_type,mode = mode,confidence = 'high').values,
        color = 'w',
        alpha = 1,
    )

plot_extreme_line(extrc,ax = ax1,extr_type='pos')
plot_extreme_line(extrc,ax = ax2,extr_type='neg')

# no grid line
# set the axis
ax1.spines["right"].set_visible(False)
ax1.spines["top"].set_visible(False)

ax1.format(
    xlabelpad=0.8, 
    xtickminor=False,
    ytickminor=False,
    xrotation=45,
    ylabel = 'ocurrence',
    xlabel = '',
    grid = False,
    ylim = (140,350),
    xlim = (ax1.get_xlim()[0] - 2000,ax1.get_xlim()[-1])
)

# same for ax2
ax2.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)

ax2.format(
    xlabelpad=0.8, 
    xtickminor=False,
    ytickminor=False,
    xrotation=45,
    ylabel = 'ocurrence',
    xlabel = '',
    grid = False,
    ylim = (140,350),
    xlim = (ax2.get_xlim()[0] - 2000,ax2.get_xlim()[-1])
)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/MPI_GE_extrc.pdf",
    dpi=300,
)


# %%
############################
# line plots for all the models

def read_extrc(model):
    """read extreme counts"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_count/"
    filename = "plev_50000_decade_mpi_first_JJA_extre_counts.nc"
    ds = xr.open_dataset(odir + filename).pc
    return ds
EXTRCs = {}
for model in ["MPI_GE", "MPI_GE_onepct", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]:
    EXTRCs[model] = read_extrc(model)


# %%
fig2 = pplt.figure(figsize=(180 / 25.4, 150 / 25.4),sharex=False,sharey=False)
fig2.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
    
)
models_legend = [
    "MPI-GE (100)",
    "MPI_GE_onepct (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    ]

gs = pplt.GridSpec(
    ncols=3,
    nrows=1,
    wspace=2,
)

ax3 = fig2.add_subplot(gs[0, 0])
ax4 = fig2.add_subplot(gs[0, 1])
ax_legend = fig2.add_subplot(gs[0, 2])

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

# put the legend in the last subplot
ax_legend.legend(
    lines_pos,
    models_legend,
    loc="c",
    ncols=1,
    frameon=False,
    pad=4,
    labelspacing = 5
)

ax_legend.spines["right"].set_visible(False)
ax_legend.spines["top"].set_visible(False)
ax_legend.spines["left"].set_visible(False)
ax_legend.spines["bottom"].set_visible(False)

ax_legend.format(
    xlabel="",
    ylabel="",
    xticklabels = [],
    yticklabels = [],
    grid = False,
    xlim = (0,1),
    ylim = (0,1),
)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/GFDL_extrc.pdf",
    dpi=300,
    bbox_inches="tight",
)


# %%

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

COMPOSITEs = {}
for model in ["MPI_GE","CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]:
    COMPOSITEs[model] = read_composite(model, "ts",reduction = 'mean')
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
params = {
    "figure.facecolor": "black",
    "axes.facecolor": "white",
    "ytick.color": "w",
    "xtick.color": "w",
    "axes.labelcolor": "w",
    "axes.edgecolor": "w",
    "tick.labelcolor": "w",
    "text.color": "w",}
pplt.rc.update(params)

axes.format(
    latlines=20,
    lonlines=30,
    color = 'grey7',
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
    leftlabels=["first10", "last10", "last10 - first10"],
    toplabels=models_legend,
    toplabels_kw = {"fontsize": 7,"color":"w"},
    leftlabels_kw = {"fontsize": 7,"color":"w"},
)

axes,maps = composite_plot.plot_composite_single_ext(COMPOSITEs, models, axes)
fig3.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K",width = 0.1,shrink=1)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/composite_pos.pdf",
    dpi=300,
    bbox_inches="tight",
)

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
params = {
    "figure.facecolor": "black",
    "axes.facecolor": "white",
    "ytick.color": "w",
    "xtick.color": "w",
    "axes.labelcolor": "w",
    "axes.edgecolor": "w",
    "tick.labelcolor": "w",
    "text.color": "w",}
pplt.rc.update(params)

axes.format(
    latlines=20,
    lonlines=30,
    color = 'grey7',
    coast=True,
    coastlinewidth=0.3,
    coastcolor="charcoal",
    leftlabels=["first10", "last10", "last10 - first10"],
    toplabels=models_legend,
    toplabels_kw = {"fontsize": 7,"color":"w"},
    leftlabels_kw = {"fontsize": 7,"color":"w"},
)

axes,maps = composite_plot.plot_composite_single_ext(COMPOSITEs, models, axes,extr_type='neg')
fig3.colorbar(maps[0], loc="b", pad=1, title=f"tsurf / K",width = 0.1,shrink=1)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/imprs_retreat/composite_neg.pdf",
    dpi=300,
    bbox_inches="tight",
)

# %%
