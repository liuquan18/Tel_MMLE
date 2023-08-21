#%%
import matplotlib.pyplot as plt
import matplotlib as mpl
import proplot as pplt
import xarray as xr
import cartopy.crs as ccrs

#%%
import src.plots.statistical_overview as stat_overview

#%%
# plot the eof spatial pattern and the positive and negative center boxes
# read the data
def read_eof_data(model,decade = 0):
    eof_name = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result.nc"
    eof_re = xr.open_dataset(eof_name)
    eof = eof_re.eof.isel(decade = decade,mode = 0)
    fra = eof_re.fra.isel(decade = decade,mode = 0)
    return eof,fra


def plot_box(lat,ulat,llon,rlon,ax,linestyle = "solid"):
    ax.plot(
        [llon,rlon,rlon,llon,llon],
        [ulat,ulat,lat,lat,ulat],
        color = "black",
        linestyle = linestyle,
        transform = ccrs.PlateCarree(),

    )
    return ax
#%%

eofs = {}
fras = {}
models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]

for model in models:
    eof,fra = read_eof_data(model)
    eofs[model] = eof
    fras[model] = fra
#%%

fig1 = pplt.figure(figsize=(180 / 25.4, 150 / 25.4),sharex=False,sharey=False)
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
    ax,fmap,_ = stat_overview.spatial_pattern_plot(
        ax,eofs[model],fras[model]
    )
    frac = ax.get_title()
    ax.set_title(f"{model} ({frac})")

    # plot the box 
    # positive box, solid line, lat,lat,lon,lon: 45, 60, -25, 5
    # negative box, dashed line, lat,lat,lon,lon: 60, 75, -70, -40

    ax = plot_box(45,60,-30,0,ax,'solid')
    ax = plot_box(60,75,-75,-50,ax,'dashed')




# %%
# read the data
models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
models_legend = [
    "MPI-GE_onepct (100)",
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
]

vars_pos = {}
vars_neg = {}

for model in models:
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

    vars_pos[model] = pos_var
    vars_neg[model] = neg_var

# %%
# plot the variability, pos as solid line, neg as dashed line
fig, axes = plt.subplots(figsize=(8, 5), ncols=2, nrows=1, sharey=True)
colors_model = ["red", "C1", "tab:purple", "tab:blue", "tab:green", "C4"]


pos_lines = []
neg_lines = []

for i, model in enumerate(models):
    pos_line = vars_pos[model].plot.line(
        x = 'time', color=colors_model[i], label=models_legend[i],
        ax=axes[0]
    )   
    neg_line = vars_neg[model].plot(
        x = 'time', color=colors_model[i], label=models_legend[i],
        ax=axes[1]
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