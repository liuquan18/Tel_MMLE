# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import proplot as pplt

# %%
from src.plots.statistical_overview import spatial_pattern_plot
from src.plots.statistical_overview import index_distribution_plot
from src.plots.extreme_plot import extrc_time_line_single
from src.plots.extreme_plot import extrc_slope_line
import seaborn as sns

# %%
from matplotlib.lines import Line2D


# %%
# read spatial patterns
def read_spatial_pattern(model, period):
    spp_path = (
        "/work/mh0033/m300883/Tel_MMLE/data/"
        + model
        + "/EOF_result/"
        + period
        + "_pattern_projected.nc"
    )
    spp = xr.open_dataset(spp_path)
    spp = spp.__xarray_dataarray_variable__
    # reaname the variable to 'NAO'
    spp = spp.rename("NAO")
    return spp


# read all spatial patterns
def read_spatial_pattern_all(model):
    eof_path = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_all_first_JJA_eof_result.nc"
    eof = xr.open_dataset(eof_path)
    pattern = eof.eof.sel(mode="NAO").isel(decade=0)
    return pattern


# %%
def read_fra(model, period):
    fra_path = (
        "/work/mh0033/m300883/Tel_MMLE/data/"
        + model
        + "/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result.nc"
    )
    fra = xr.open_dataset(fra_path)
    fra = fra.fra.sel(mode="NAO")
    if period == "first":
        fra = fra.isel(decade=0)
    elif period == "last":
        fra = fra.isel(decade=-1)
    return fra


def read_fra_all(model):
    eof_path = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_all_first_JJA_eof_result.nc"
    eof = xr.open_dataset(eof_path)
    fra = eof.fra.sel(mode="NAO").isel(decade=0)
    return fra


# %%
def read_pc(model, period, fixed_pattern="decade_mpi"):
    pc_path = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_{fixed_pattern}_first_JJA_eof_result.nc"
    pc = xr.open_dataset(pc_path)
    pc = pc.pc.sel(mode="NAO")
    times = pc.time
    years = np.unique(times.dt.year)
    if period == "first":
        pc = pc.sel(time=pc["time.year"].isin(years[:10]))
    elif period == "last":
        pc = pc.sel(time=pc["time.year"].isin(years[-10:]))
    return pc


# %%
def read_extrc_fixed(model):
    extrc_path = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/plev_50000_all_first_JJA_extre_counts.nc"
    extrc = xr.open_dataset(extrc_path)
    extrc = extrc.pc

    ens_sizes = {
        "CanESM2": 50,
        "CESM1_CAM5": 40,
        "MK36": 30,
        "GFDL_CM3": 20,
        "MPI_GE_onepct": 100,
        "MPI_GE": 100,
    }

    # divide by ensemble size
    extrc = extrc / ens_sizes[model]

    return extrc


# %%
def read_slope_fixed(model):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_slope/"
    filename = f"plev_50000_all_first_JJA_extre_slope_time.nc"
    ds = xr.open_dataset(odir + filename).pc
    return ds


# %%
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


# %%
def read_ens_var(model):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/ens_variability/{model}_ens_var_10y.nc"
    ds = xr.open_dataset(odir)
    ens_var = ds.zg
    return ens_var


# %%
First_patterns = {}
Last_patterns = {}

All_patterns = {}  # the patterns from the whole length of the simulation
All_fras = (
    {}
)  # the fraction of explained variance from the whole length of the simulation

First_fras = {}
Last_fras = {}

First_pcs_decade = {}
Last_pcs_decade = {}

First_pcs_all = {}
Last_pcs_all = {}

EXTRCs = {}
EXTRCs_slope = {}


# SMILEs
for model in ["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "MPI_GE_onepct", "MPI_GE"]:
    First_patterns[model] = read_spatial_pattern(model, "first")
    Last_patterns[model] = read_spatial_pattern(model, "last")

    All_patterns[model] = read_spatial_pattern_all(model)
    All_fras[model] = read_fra_all(model)

    First_fras[model] = read_fra(model, "first")
    Last_fras[model] = read_fra(model, "last")

    First_pcs_decade[model] = read_pc(model, "first", "decade_mpi")
    Last_pcs_decade[model] = read_pc(model, "last", "decade_mpi")

    First_pcs_all[model] = read_pc(model, "first", "all")
    Last_pcs_all[model] = read_pc(model, "last", "all")

    EXTRCs[model] = read_extrc_fixed(model)
    EXTRCs_slope[model] = read_slope_fixed(model)
# %%
EXTRCs_rates = SMILE_rate_increase(EXTRCs_slope, EXTRCs)

# %%
# CR20_allens
CR20_first_pattern = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_first_40_pattern_projected.nc"
)
CR20_last_pattern = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_last_40_pattern_projected.nc"
)
# %%
First_patterns["CR20_allens"] = CR20_first_pattern.NAO
Last_patterns["CR20_allens"] = CR20_last_pattern.NAO
# %%

First_fras["CR20_allens"] = (
    xr.open_dataset(
        "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_first_40_eof_std.nc"
    )
    .fra.sel(mode="NAO")
    .values[0]
)
Last_fras["CR20_allens"] = (
    xr.open_dataset(
        "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_last_40_eof_std.nc"
    )
    .fra.sel(mode="NAO")
    .values[0]
)

# %%
First_pcs_decade["CR20_allens"] = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_first_40_eof_std.nc"
).pc.sel(mode="NAO")
Last_pcs_decade["CR20_allens"] = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_last_40_eof_std.nc"
).pc.sel(mode="NAO")
# %%
First_pcs_all["CR20_allens"] = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/first_40_eof_std.nc"
).pc.sel(mode="NAO")
Last_pcs_all["CR20_allens"] = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/last_40_eof_std.nc"
).pc.sel(mode="NAO")

# %%
################# Fig.S1 spatial patterns #################
fig1, axes = pplt.subplots(
    ncols=3,
    nrows=2,
    proj="ortho",
    proj_kw={"lon_0": -20, "lat_0": 60},
    abc=True,
    figsize=(180 / 25.4, 150 / 25.4),
)
models_legend = [
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    "MPI_GE_onepct (100)",
    "CR20_ensemble (40)",
]

models = ["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "MPI_GE_onepct", "CR20_allens"]

# flat the axes so that we can iterate over them
axes = np.array(axes).ravel()

for i, ax in enumerate(axes):
    model = models[i]
    spatial_ax, fmap, lmap = spatial_pattern_plot(
        ax,
        First_patterns[model],
        First_fras[model],
        Last_patterns[model],
        Last_fras[model],
        levels=np.arange(-30, 31, 5),
    )
fig1.colorbar(fmap, loc="b", label="500hPa geopotential height (m)")

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/figS1_spatial_pattern.pdf",
)

# %%
################# Fig.S2 index distribution #################

fig2, axes = pplt.subplots(
    ncols=3,
    nrows=2,
    abc=True,
    share=True,
    figsize=(180 / 25.4, 150 / 25.4),
)


models_legend = [
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    "MPI_GE_onepct (100)",
    "CR20_ensemble (40)",
]

models = ["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "MPI_GE_onepct", "CR20_allens"]

# flat the axes so that we can iterate over them
axes = np.array(axes).ravel()

for i, ax in enumerate(axes):
    model = models[i]
    index_distribution_plot(
        ax,
        First_pcs_decade[model],
        Last_pcs_decade[model],
    )

    first_std = First_pcs_decade[model].std().values
    last_std = Last_pcs_decade[model].std().values

    ax.set_title(
        "({first_std:0.2f} -> {last_std:0.2f})".format(
            first_std=first_std, last_std=last_std
        )
    )

    ax.set_ylabel("Probability density")
    ax.set_xlabel("NAO index")

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/figS2_index_distribution.pdf",
)

# %%
################# Fig.S3 spatial patterns from the whole periods #################
fig3, axes = pplt.subplots(
    ncols=3,
    nrows=2,
    proj="ortho",
    proj_kw={"lon_0": -20, "lat_0": 60},
    abc=True,
    figsize=(180 / 25.4, 150 / 25.4),
)
models_legend = [
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    "MPI_GE_onepct (100)",
    "MPI_GE (100)",
]

models = ["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "MPI_GE_onepct", "MPI_GE"]

# flat the axes so that we can iterate over them
axes = np.array(axes).ravel()

for i, ax in enumerate(axes):
    model = models[i]
    spatial_ax, fmap, lmap = spatial_pattern_plot(
        ax,
        All_patterns[model],
        All_fras[model],
        levels=np.arange(-30, 31, 5),
    )
fig3.colorbar(fmap, loc="b", label="500hPa geopotential height (m)")
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/figS3_spatial_pattern_all.pdf",
)

# %%
################# Fig.S4 index distribution for 'all' pattern #################

fig4, axes = pplt.subplots(
    ncols=3,
    nrows=2,
    abc=True,
    share=True,
    figsize=(180 / 25.4, 150 / 25.4),
)


models_legend = [
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
    "MPI_GE_onepct (100)",
    "MPI_GE (100)",
]

models = ["CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "MPI_GE_onepct", "MPI_GE"]

# flat the axes so that we can iterate over them
axes = np.array(axes).ravel()

for i, ax in enumerate(axes):
    model = models[i]
    index_distribution_plot(
        ax,
        First_pcs_all[model],
        Last_pcs_all[model],
    )

    first_std = First_pcs_all[model].std().values
    last_std = Last_pcs_all[model].std().values

    ax.set_title(
        "({first_std:0.2f} -> {last_std:0.2f})".format(
            first_std=first_std, last_std=last_std
        )
    )

    ax.set_ylabel("Probability density")
    ax.set_xlabel("NAO index")

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/figS4_index_distribution_all.pdf",
)


# %%
################# Fig.S5 extreme counts #################
fig5, axes = pplt.subplots(
    ncols=2, nrows=2, abc=True, share=False, figsize=(150 / 25.4, 150 / 25.4)
)
extrc_time_line_single(
    EXTRCs,
    extr_type="pos",
    ax=axes[0, 0],
    ylim=(1, 3.5),
    ci=True,
)

extrc_time_line_single(
    EXTRCs,
    extr_type="neg",
    ax=axes[0, 1],
    ylim=(1, 3.5),
    ci=True,
)

extrc_slope_line(
    EXTRCs_rates,
    extr_type="pos",
    ax=axes[1, 0],
    models=["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"],
)

extrc_slope_line(
    EXTRCs_rates,
    extr_type="neg",
    ax=axes[1, 1],
    models=["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"],
)

ax1, ax2, ax3, ax4 = np.array(axes).ravel()

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
    ylim=(-2.4, 8),
    xlim=(16, 72),
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

ax3.spines["bottom"].set_position(("data", 0))
ax4.spines["bottom"].set_position(("data", 0))

ax4.format(
    ytickminor=False,
    grid=False,
    ylabel="",
    ylim=(-2.4, 8),
    xtickminor=False,
    facecolor="none",
    xlim=(16, 72),
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
]


fig5.legend(
    legend_lines,
    models_legend,
    loc="b",
    frameon=False,
)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/figS5_extreme_counts_all.pdf",
)

# %%
#################### Fig.S6 20CR v.s MPI_GE ####################

# read the pc of MPI_GE with all pattern
MPI_GE_all = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/plev_50000_all_first_JJA_eof_result.nc"
)
MPI_GE_pc = MPI_GE_all.pc.sel(mode="NAO", time=slice("1850", "2015"))
MPI_GE_pc.name = "NAO index"

# read the pc of 20CR
CR20_all = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/all_eof.nc"
)
CR20_pc = CR20_all.pc.sel(mode="NAO", time=slice("1850", "2015"))
CR20_pc.name = "NAO index"

# %%
# read the mean zg at the south center
MPI_GE_zg_south = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/ens_variability/MPI_GE_zg_south_center.nc"
)
MPI_GE_zg_south = MPI_GE_zg_south.var156
MPI_GE_zg_south = MPI_GE_zg_south.sel(time=slice("1850", "2015"))
MPI_GE_zg_south.name = "500hPa geopotential height (m)"

# %%
CR20_zg_south = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/ens_variability/CR20_zg_south_center.nc"
)
CR20_zg_south = CR20_zg_south.HGT
CR20_zg_south.name = "500hPa geopotential height (m)"

# %%
# read the ens_var
MPI_GE_ens_var = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/ens_variability/MPI_GE_ens_var_10y.nc"
)
MPI_GE_ens_var = MPI_GE_ens_var.zg
MPI_GE_ens_var.name = "Ensemble standard deviation (m)"
# %%
CR20_ens_var = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/ens_variability/CR20_ens_var_10y.nc"
)
CR20_ens_var = CR20_ens_var.zg
CR20_ens_var.name = "Ensemble standard deviation (m)"
# %%
CR20_ens_var_40 = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/ens_variability/CR20_ens_var_40y.nc"
)
CR20_ens_var_40 = CR20_ens_var_40.__xarray_dataarray_variable__
CR20_ens_var_40.name = "Ensemble standard deviation (m)"
# %%
MPI_GE_ens_var_40 = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/ens_variability/MPI_GE_ens_var_40y.nc"
)
MPI_GE_ens_var_40 = MPI_GE_ens_var_40.__xarray_dataarray_variable__
MPI_GE_ens_var_40.name = "Ensemble standard deviation (m)"
# %%
fig6, axes = pplt.subplots(
    ncols=2,
    nrows=3,
    abc=True,
    sharex=True,
    sharey=False,
    figsize=(150 / 25.4, 180 / 25.4),
)

MPI_GE_zg_south.plot.line(x="time", ax=axes[0, 0], lw=1, add_legend=False)
CR20_zg_south.plot.line(x="time", ax=axes[0, 1], lw=1, add_legend=False)

MPI_GE_ens_var.plot.line(x="time", ax=axes[1, 0], color="k", lw=1.0, add_legend=False)
MPI_GE_ens_var_40.plot.line(
    x="time", ax=axes[1, 0], color="k", lw=1.0, add_legend=False, linestyle="--"
)

CR20_ens_var.plot.line(
    x="time",
    ax=axes[1, 1],
    color="k",
    lw=1.0,
    add_legend=False,
)
CR20_ens_var_40.plot.line(
    x="time", ax=axes[1, 1], color="k", lw=1.0, add_legend=False, linestyle="--"
)

MPI_GE_pc.plot.line(x="time", ax=axes[2, 0], lw=1.0, add_legend=False)
CR20_pc.plot.line(
    x="time",
    ax=axes[2, 1],
    lw=1.0,
    add_legend=False,
)


axes.format(
    grid=False,
    ytickminor=False,
    xtickminor=False,
    facecolor="none",
)

axes[:, 1].set_ylabel("")
axes[0, :].set_ylim(5580, 5851)
axes[1, :].set_ylim(20.9, 30.9)
axes[2, :].set_ylim(-4, 4.3)

for ax in axes:
    ax.set_title("")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/figS6_20CR_vs_MPI_GE.pdf",
)


# %%
###################### direct NAO response ######################
def read_index(model):
    file_path = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_all_first_JJA_projected_index.nc"
    index = xr.open_dataset(file_path).pc.sel(mode="NAO")
    return index


# %%
def cal_slope(index):
    index_notime = index.assign_coords(time=np.arange(len(index["time"])))
    # calculate the linear slope of mpi_index along time dimension
    res = index_notime.polyfit(dim="time", deg=1)
    slope = res.polyfit_coefficients.sel(degree=1)
    return slope


# %%
def convert2pd(model):
    index = read_index(model)
    slope = cal_slope(index)
    # every decade
    slope = slope * 10
    return slope.to_dataframe().reset_index()[["polyfit_coefficients"]]


# %%

mpi_slope = convert2pd("MPI_GE")
canesm2_slope = convert2pd("CanESM2")
cesm1_cam5_slope = convert2pd("CESM1_CAM5")
gfdl_cm3_slope = convert2pd("GFDL_CM3")
mk36_slope = convert2pd("MK36")

# %%
# combine all dataframes into one
combined_df = pd.concat(
    [mpi_slope, canesm2_slope, cesm1_cam5_slope, gfdl_cm3_slope, mk36_slope], axis=1
)

# %%
combined_df.columns = ["MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
# %%

# make box plot
plt.figure(figsize=(100 / 25.4, 100 / 25.4))
sns.boxplot(
    data=combined_df, palette=["C1", "tab:purple", "tab:blue", "tab:green", "yellow"]
)
plt.xticks(
    range(5),
    ["MPI_GE(100)", "CanESM2(50)", "CESM1_CAM5(40)", "GFDL_CM3(30)", "MK36(20)"],
)
plt.ylabel("Slope / 10 yr")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/direct_response.pdf"
)

# %%
################ zonal wind v.s blocking ###################


# %%
def read_blocking(model):
    fpath = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/GB_index_JJA"
    blocking = xr.open_mfdataset(
        f"{fpath}/trend*.nc", combine="nested", concat_dim="ens"
    )
    blocking = blocking.squeeze()
    try:
        blocking = blocking["zg"].to_dataframe().reset_index()[["zg"]]
    except KeyError:
        try:
            blocking = blocking["gph"].to_dataframe().reset_index()[["gph"]]
        except KeyError:
            blocking = blocking["var156"].to_dataframe().reset_index()[["var156"]]
    blocking.columns = ["blocking"]
    return blocking.dropna()


# %%
def read_zonal(model):
    fpath = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/zonal_wind_JJA"
    zonal = xr.open_mfdataset(f"{fpath}/trend*.nc", combine="nested", concat_dim="ens")
    zonal = zonal.squeeze()
    try:
        zonal = zonal["var131"].to_dataframe().reset_index()[["var131"]]
    except KeyError:
        zonal = zonal["ua"].to_dataframe().reset_index()[["ua"]]
    zonal.columns = ["zonal_wind"]
    return zonal.dropna()


# %%
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]  #'CESM1_CAM5''MPI_GE',
# %%
# concant all models into one dataframe
diagnose = {}
for model in models:
    diagnose[model] = pd.concat([read_blocking(model), read_zonal(model)], axis=1)


# %%
# sns scatter plot, x-axis blocking, y-axis zonal wind
fig, axes = plt.subplots(
    2, 3, figsize=(180 / 25.4, 180 / 25.4), sharex=True, sharey=True
)
axes = axes.flatten()
models_plot = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
colors = ["C1", "tab:purple", "tab:blue", "tab:green", "yellow"]
for i, model in enumerate(models_plot):
    sns.scatterplot(
        data=diagnose[model], x="blocking", y="zonal_wind", ax=axes[i], color=colors[i]
    )
    axes[i].set_title(model)
    axes[i].set_xlabel("")
    axes[i].set_ylabel("")

X = [diagnose[model].blocking.mean() for model in models_plot]
Y = [diagnose[model].zonal_wind.mean() for model in models_plot]
Xerr = [diagnose[model].blocking.std() for model in models_plot]
Yerr = [diagnose[model].zonal_wind.std() for model in models_plot]
for i, model in enumerate(models_plot):
    axes[-1].errorbar(
        X[i],
        Y[i],
        xerr=Xerr[i],
        yerr=Yerr[i],
        fmt="o",
        color=colors[i],
        label=model,
        marker="o",
    )
axes[-1].set_title("All models")
for ax in axes:
    ax.axhline(0, color="black", lw=0.5)
    ax.axvline(0, color="black", lw=0.5)
# add the y-lable "zonal wind (200 hPa)" for the fig
fig.text(0.04, 0.5, "Zonal wind (200 hPa)", va="center", rotation="vertical", fontsize = 10)
# add the x-lable "Greenland blocking index (500 hPa)" for the fig
fig.text(0.5, 0.04, "Greenland blocking index (500 hPa)", ha="center", fontsize = 10)

legend_lines = [
    Line2D([0], [0], color=colors[0], lw=1.5, marker="o"),
    Line2D([0], [0], color=colors[1], lw=1.5, marker="o"),
    Line2D([0], [0], color=colors[2], lw=1.5, marker="o"),
    Line2D([0], [0], color=colors[3], lw=1.5, marker="o"),
    Line2D([0], [0], color=colors[4], lw=1.5, marker="o"),
]


fig.legend(
    legend_lines,
    models_plot,
    loc="lower center",
    frameon=False,
    ncol=5,
)
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/blocking_zonal_wind.pdf"
)
# %%
