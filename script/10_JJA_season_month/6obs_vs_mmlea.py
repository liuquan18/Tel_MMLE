# %%
import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import proplot as pplt
import matplotlib.ticker as ticker

# %%
import src.plots.statistical_overview as sov
import importlib

importlib.reload(sov)


# %%
def read_eof_result(model, time=slice("1940", "2022")):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/ERA5/EOF_result/"
    filename = f"plev_50000_1940_2022_{model}_all.nc"
    ds = xr.open_dataset(odir + filename)
    NAO = ds.sel(mode="NAO").isel(decade=0)
    NAO = NAO.sortby("time").sel(time=time)
    return NAO

# %%
def obs_mmlea_spatial_pattern(EOFs):
    fig = pplt.figure(space=0, refwidth="20em", wspace=3, hspace=3)
    gs = pplt.GridSpec(ncols=3,nrows=2)

    for i, model in enumerate(EOFs.keys()):
        spatial_ax, fmap, lmap = sov.spatial_pattern_plot(
            fig,
            gs[divmod(i, 3)],
            EOFs[model].eof,
            EOFs[model].fra,
            levels=np.arange(-40, 41, 5),
            title=models_legend[i],
        )
    fig.colorbar(
        fmap,
        loc="r",
    )
    plt.savefig(
        "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/obs_vs_MMLEA/EOF_spatial_pattern.png",
    )


# %%
def obs_mmlea_time_series(EOFs):
    fig = pplt.figure(
        space=0, refwidth="25em", wspace=3, hspace=3, sharex=False, sharey=False
    )

    gs = pplt.GridSpec(ncols=3,nrows=2)
    for i, model in enumerate(models):
        obs = EOFs["ERA5"].pc
        mmlea = EOFs[model].pc
        # the line plot
        sov.envelop_obs_mmlea(fig, gs[divmod(i, gs.ncols)], obs, mmlea)

    # the box plot
    box_ax = fig.add_subplot(gs[1, 2])
    sov.obs_mmlea_box_plot(ax,EOFs)
    plt.savefig(
        "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/obs_vs_MMLEA/EOF_time_series.png",
    )


# %%
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
models_legend = [
    "ERA5 (obs)",
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
]

EOFs = {}
EOFs["ERA5"] = read_eof_result("ERA5", time=slice("1960", "2022"))
for model in models:
    EOFs[model] = read_eof_result(model, time=slice("1960", "2022"))
# %%
obs_mmlea_spatial_pattern(EOFs)
# %%
obs_mmlea_time_series(EOFs)
# %%
