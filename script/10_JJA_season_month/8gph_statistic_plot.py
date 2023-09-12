# %%
import matplotlib.pyplot as plt
import matplotlib as mpl
import proplot as pplt
import xarray as xr
import cartopy.crs as ccrs
import numpy as np

# %%
import src.plots.statistical_overview as stat_overview


# %%
# plot the eof spatial pattern and the positive and negative center boxes
# read the data
def read_eof_data(model, decade=0):
    eof_name = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result.nc"
    eof_re = xr.open_dataset(eof_name)
    eof = eof_re.eof.isel(mode=0)
    fra = eof_re.fra.isel(mode=0)
    return eof, fra


def read_box_stats(model):
    pos_name = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/box_based/pos_var.nc"
    neg_name = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/box_based/neg_var.nc"

    pos_var = xr.open_dataset(pos_name)
    neg_var = xr.open_dataset(neg_name)

    try:
        pos_var = pos_var.var156
        neg_var = neg_var.var156
    except AttributeError:
        pos_var = pos_var.zg
        neg_var = neg_var.zg
    return pos_var, neg_var


def read_slope_stat(model,var='end_std'):
    slope_name = (
        f"/work/mh0033/m300883/Tel_MMLE/data/{model}/box_based/slope_of_{var}.nc"
    )
    slope = xr.open_dataset(slope_name)
    return slope



# %%
# read data for all the models
eofs = {}
fras = {}

vars_pos = {}  # for the variability chagne over boxes
vars_neg = {}

slopes_ens_std = {}  # the slope of the ensemble std over the whole North Atlantic sector
slopes_ens_mean = {}
slopes_ens_extrc = {}

models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]

# %%
for model in models:
    eof, fra = read_eof_data(model)
    eofs[model] = eof
    fras[model] = fra

    pos_var, neg_var = read_box_stats(model)
    vars_pos[model] = pos_var
    vars_neg[model] = neg_var

    slope_std = read_slope_stat(model, var='ens_std')
    slopes_ens_std[model] = slope_std

    slope_mean = read_slope_stat(model, var='ens_mean')
    slopes_ens_mean[model] = slope_mean


    slope_extrc = read_slope_stat(model, var='gph_extrc')
    slopes_ens_extrc[model] = slope_extrc


# %%
# %%
models_legend = [
    "MPI-GE_onepct (100)",
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
]


def plot_box_outline(lat, ulat, llon, rlon, ax, linestyle="solid"):
    ax.plot(
        [llon, rlon, rlon, llon, llon],
        [ulat, ulat, lat, lat, ulat],
        color="black",
        linestyle=linestyle,
        transform=ccrs.PlateCarree(),
    )
    return ax


def plot_slope_std_singleModel(ax, slope,pvalue, sig = False,levels = np.arange(-0.8,0.9,0.2)):
    map = slope.plot(
        ax=ax,
        color="black",
        linestyle="solid",
        transform=ccrs.PlateCarree(),
        add_colorbar=False,
        levels = levels,
        extend = 'both'
    )

    if sig:
        sig = slope.plot.contourf(
            ax=ax,
            levels=[0, 0.05],
            colors="none",
            hatches=["", "///"],
            transform=ccrs.PlateCarree(),
            add_colorbar=False,
        )
    
    return ax, map

def plot_slope_mean_singleModel(ax, slope, sig = False):
    map = ax.contour(
        slope.lon,
        slope.lat,
        slope.slope,
        levels = np.arange(4,14,0.5),
        colors = 'grey8',
        nozero = True,
        add_colorbar = False,
        lw = 0.8,
        labels = True,
        transform = ccrs.PlateCarree(),
    )

    if sig:
        sig = slope.pvalue.plot.contourf(
            ax=ax,
            levels=[0, 0.05],
            colors="none",
            hatches=["", "///"],
            transform=ccrs.PlateCarree(),
            add_colorbar=False,
        )
    
    return ax, map


# %%

fig1 = pplt.figure(figsize=(180 / 25.4, 150 / 25.4), sharex=False, sharey=False)
fig1.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

axes = fig1.subplots(
    ncols=3,
    nrows=2,
    proj="ortho",
    proj_kw={"lon_0": -20, "lat_0": 60},
)

for i, ax in enumerate(axes):
    model = models[i]
    ax, fmap, _ = stat_overview.spatial_pattern_plot(ax, eofs[model].isel(decade = 0), 
                                                     fras[model].isel(decade = 0),
                                                     eofs[model].isel(decade = -1),
                                                     fras[model].isel(decade = -1),)
    frac = ax.get_title()
    ax.set_title(f"{model} ({frac})")

    # plot the box
    # positive box, solid line, lat,lat,lon,lon: 45, 60, -25, 5
    # negative box, dashed line, lat,lat,lon,lon: 60, 75, -70, -40

    ax = plot_box_outline(45, 60, -30, 0, ax, "solid")
    ax = plot_box_outline(60, 75, -75, -50, ax, "dashed")

fig1.colorbar(fmap, loc="b", label="NAO standard", shrink=0.8)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/box_loc.png"
)

#%%
# plot the variability of the positive and negative center over time
fig2, axes = plt.subplots(figsize=(8, 5), ncols=2, nrows=1, sharey=True)
colors_model = ["red", "C1", "tab:purple", "tab:blue", "tab:green", "C4"]


pos_lines = []
neg_lines = []

for i, model in enumerate(models):
    pos_line = vars_pos[model].plot.line(
        x="time", color=colors_model[i], label=models_legend[i], ax=axes[0]
    )
    neg_line = vars_neg[model].plot(
        x="time", color=colors_model[i], label=models_legend[i], ax=axes[1]
    )
    pos_lines.append(pos_line)
    neg_lines.append(neg_line)

axes[0].set_title("Positive center")
axes[1].set_title("Negative center")

axes[1].legend(loc="best")

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/variability_pos_neg_center.png"
)

# %%
# slope of the variability along ensemble dimension over time
fig3 = pplt.figure(figsize=(180 / 25.4, 150 / 25.4), sharex=False, sharey=False)
fig3.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

axes = fig3.subplots(
    ncols=3,
    nrows=2,
    proj="ortho",
    proj_kw={"lon_0": -20, "lat_0": 60},
)

for i, ax in enumerate(axes):
    model = models[i]
    ax = plot_box_outline(45, 60, -30, 0, ax, "solid")
    ax = plot_box_outline(60, 75, -75, -50, ax, "dashed")
    ax, map = plot_slope_std_singleModel(ax=axes[i], slope=slopes_ens_std[model].slope, pvalue=slopes_ens_std[model].pvalue, sig = True)
    ax.format(
        lonlines=20,
        latlines=30,
        coast=True,
        coastlinewidth=0.5,
        coastcolor="charcoal",
        title=f"{models_legend[i]}",
    )

fig3.colorbar(map, 
              orientation="horizontal", 
              shrink=0.5, 
              loc="b",
              title = 'slope of the variability along ensemble dimension every ten years')
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/slope_ens_std.png"
)


# %%
# Fig 4: slope of the std and mean along ensemble dimension over time
fig4 = pplt.figure(figsize=(180 / 25.4, 150 / 25.4), sharex=False, sharey=False)
fig4.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

axes = fig4.subplots(
    ncols=3,
    nrows=2,
    proj="ortho",
    proj_kw={"lon_0": -20, "lat_0": 60},
)

for i, ax in enumerate(axes):
    model = models[i]

    try:
        ax, map_std = plot_slope_std_singleModel(ax=axes[i], slope=slopes_ens_std[model].slope, pvalue=slopes_ens_std[model].pvalue)
    except:
        pass
    try:
        ax, map_mean = plot_slope_mean_singleModel(ax=axes[i], slope=slopes_ens_mean[model])
    except:
        pass
    ax.format(
        lonlines=20,
        latlines=30,
        coast=True,
        coastlinewidth=0.5,
        coastcolor="charcoal",
        title=f"{models_legend[i]}",
    )

fig4.colorbar(map_std, 
              orientation="horizontal", 
              shrink=0.5, 
              loc="b",
              title = 'slope of the variability along ensemble dimension every ten years')

fig4.colorbar(map_mean,
                orientation="horizontal", 
                shrink=0.5, 
                loc="b",
                title = 'slope of the mean along ensemble dimension every ten years')

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/slope_ens_mean.png"
)


# %%
# Fig 5: slope of the pos extreme events along ensemble dimension over time
fig5 = pplt.figure(figsize=(180 / 25.4, 150 / 25.4), sharex=False, sharey=False)
fig5.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

axes = fig5.subplots(
    ncols=3,
    nrows=2,
    proj="ortho",
    proj_kw={"lon_0": -20, "lat_0": 60},
)

for i, ax in enumerate(axes):
    model = models[i]
    ax = plot_box_outline(45, 60, -30, 0, ax, "solid")
    ax = plot_box_outline(60, 75, -75, -50, ax, "dashed")
    ax, map = plot_slope_std_singleModel(ax=axes[i], 
                                        slope=slopes_ens_extrc[model].slope_pos,
                                        pvalue=slopes_ens_extrc[model].pvalue_pos,
                                        levels = np.arange(-4,5,1),)
    ax.format(
        lonlines=20,
        latlines=30,
        coast=True,
        coastlinewidth=0.5,
        coastcolor="charcoal",
        title=f"{models_legend[i]}",
    )

fig5.colorbar(map, 
              orientation="horizontal", 
              shrink=0.5, 
              loc="b",
              title = 'slope of the extreme occurence every ten years')

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/slope_pos_extrc.png"
)

# %%
# Fig 5: slope of the neg extreme events along ensemble dimension over time
fig6 = pplt.figure(figsize=(180 / 25.4, 150 / 25.4), sharex=False, sharey=False)
fig6.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
)

axes = fig6.subplots(
    ncols=3,
    nrows=2,
    proj="ortho",
    proj_kw={"lon_0": -20, "lat_0": 60},
)

for i, ax in enumerate(axes):
    model = models[i]
    ax = plot_box_outline(45, 60, -30, 0, ax, "solid")
    ax = plot_box_outline(60, 75, -75, -50, ax, "dashed")
    ax, map = plot_slope_std_singleModel(ax=axes[i], 
                                        slope=slopes_ens_extrc[model].slope_neg,
                                        pvalue=slopes_ens_extrc[model].pvalue_neg,
                                        levels = np.arange(-4,5,1),)
    ax.format(
        lonlines=20,
        latlines=30,
        coast=True,
        coastlinewidth=0.5,
        coastcolor="charcoal",
        title=f"{models_legend[i]}",
    )

fig6.colorbar(map, 
              orientation="horizontal", 
              shrink=0.5, 
              loc="b",
              title = 'slope of the extreme occurence every ten years')

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/slope_neg_extrc.png"
)

# %%
