# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from matplotlib.lines import Line2D
import matplotlib.patches as mpatches


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
models = ["MPI_GE", "CanESM2",'CESM1_CAM5', "MK36", "GFDL_CM3"]  #'CESM1_CAM5''MPI_GE',
# %%
# concant all models into one dataframe
diagnose = {}
for model in models:
    diagnose[model] = pd.concat([read_blocking(model), read_zonal(model)], axis=1)


#%%
diagnose['CESM1_CAM5']['blocking'] = read_blocking('CESM1_CAM5')['blocking']

# %%
# sns scatter plot, x-axis blocking, y-axis zonal wind
fig, axes = plt.subplots(2, 3, figsize=(180 / 25.4, 180 / 25.4), sharex=True, sharey=True)
axes = axes.flatten()
models_plot = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
colors = ["C1", "tab:purple", "tab:blue", "tab:green", "yellow"]
for i, model in enumerate(models_plot):
    sns.scatterplot(
        data=diagnose[model], x="blocking", y="zonal_wind", ax=axes[i], color=colors[i]
    )
    axes[i].set_title(model)

X = [diagnose[model].blocking.mean() for model in models_plot]
Y = [diagnose[model].zonal_wind.mean() for model in models_plot]
Xerr = [diagnose[model].blocking.std() for model in models_plot]
Yerr = [diagnose[model].zonal_wind.std() for model in models_plot]
for i, model in enumerate(models_plot):
    axes[-1].errorbar(
        X[i], Y[i], xerr=Xerr[i], yerr=Yerr[i], fmt="o", color=colors[i], label=model
    )
for ax in axes:
    ax.axhline(0, color="black", lw=0.5)
    ax.axvline(0, color="black", lw=0.5)

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
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/blocking_zonal_wind.pdf")
# %%
