#%%
import xarray as xr
import numpy as np
import src.extreme.extreme_ci as extreme
import proplot as pplt
import warnings

#%%
def extCount_tsurf_scatter(
    ext_counts, t_surf, ylim=(0, 55), xlim=(-1, 5), xlabel="temperature (K)"
):
    """
    rows: pos/neg
    cols: NAO/EA
    scatter: extreme_count v.s surface temperature
    hue: different  dataset
    """
    fig, axes = pplt.subplots(nrows=2, ncols=2, figwidth=8, span=False, share=False)

    axes.format(
        abc="a",
        abcloc="ul",
        # xlim=xlim,
        suptitle=f"extreme counts v.s surface temperature in decadal time scale",
        xlabel=xlabel,
        ylabel="extreme count",
        grid=False,
        leftlabels=["NAO", "EA"],
        toplabels=["pos", "neg"],
        xminorticks="null",
        yminorticks="null",
        # ylim=ylim,
    )

    scatter_plot(ext_counts, t_surf,  axes)


def scatter_plot(ext_counts, t_surf,  axes):
    for i, mode in enumerate(ext_counts.mode):
        for j, extr_type in enumerate(ext_counts.extr_type):
            # data preparation
            true = ext_counts.sel(extr_type=extr_type, mode=mode, confidence="true")
            low = ext_counts.sel(extr_type=extr_type, mode=mode, confidence="low")
            high = ext_counts.sel(extr_type=extr_type, mode=mode, confidence="high")

            # for the data with plev
            try:
                true = true.stack(com=("time", "plev"))
                low = low.stack(com=("time", "plev"))
                high = high.stack(com=("time", "plev"))
                t = t_surf.stack(com=("time", "plev"))

            except KeyError:
                t = t_surf

            axes[i, j].errorbar(
                x=t,
                y=true,
                yerr=[(true - low), (high - true)],
                fmt="o",
                linewidth=2,
                capsize=6,
            )

            axes[i, j].set_xlim(-1, 5)

