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
    models = ["MPI_GE_onepct", "MPI_GE"]#, "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
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
plt.rcParams['hatch.linewidth'] = 0.3
# %%
EOFs_decade = read_all_models("eof_decade")
EOFs_all = read_all_models("eof_all")
EXTRCs = read_all_models("extrc")
COMPOSITEs = read_all_models("composite")

#%%

def extrc_time_line_single(extrcs, extr_type, ax, ylim = (20, 280),mode = 'NAO',ci = False):
    models = ["MPI_GE", "MPI_GE_onepct"]#,"CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    colors_model = ["C1", "red"]#, "tab:purple", "tab:blue", "tab:green", "C4"]
    model_color = dict(zip(models, colors_model))

    lines = []
    for model in models:
        try:
            extrc = extrcs[models.index(model)]
        except KeyError:
            extrc = extrcs[model]
        line = extrc.sel(extr_type=extr_type,mode = mode,confidence = 'true').plot.line(
                ax=ax, 
                label=model_color[model],
                x = 'time',
                color = model_color[model],
                linewidth = 2,)
        
        if ci:
            # fill between the confidence interval ['low','high']
            ax.fill_between(
                extrc.time,
                extrc.sel(extr_type=extr_type,mode = mode,confidence = 'low').values,
                extrc.sel(extr_type=extr_type,mode = mode,confidence = 'high').values,
                color = model_color[model],
                alpha = 0.3,
            )
        lines.append(line)

            
    ax.format(
            ylim=ylim,
            ylabel="Extreme counts",
            xlabel="Year",
            suptitle="",
            titleloc="uc",
            ylocator=20,
            yminorlocator="null",
            grid=False,
            title = '',
        )
    return ax, lines

# %%
# Fig 2, spatial pattern chagne, index distribution change, and extreme lines
# set the fig size as 180mm x 150mm
# set the font size as 7pt
fig2 = pplt.figure(figsize=(100 / 25.4, 100 / 25.4),sharex=False,sharey=False)
fig2.format(
    abc=True,
)

gs = pplt.GridSpec(
    ncols=2,
    nrows=1,
    wspace=1,
)

ax3 = fig2.add_subplot(gs[0, 0])
ax4 = fig2.add_subplot(gs[0, 1])


ax3,lines_pos = extrc_time_line_single(
    EXTRCs,
    extr_type='pos',
    ax = ax3,
    ylim = (130,330),
    ci = True,
)

ax4,lines_neg = extrc_time_line_single(
    EXTRCs,
    extr_type='neg',
    ax = ax4,
    ylim = (130,330),
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
    title = 'positive'
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
    title = 'negative'
)

ax4.legend(
    lines_pos,
    ['MPI_GE_Hist_RCP8.5(100)','MPI_GE_onepct(100)'],
    loc="b",
    ncols=2,
    frameon=False,
    bbox_to_anchor=(-0.1, -0.2),
    pad = -0.2,
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
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/composite_tsurf_pos.png",
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
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Story_line_nature_climate_change/composite_tsurf_neg.png",
    dpi=300,
    bbox_inches="tight",
)

# %%
