#%%
import xarray as xr
import numpy as np
import src.extreme.extreme_ci as extreme
import proplot as pplt

# reimport extreme
import importlib
importlib.reload(extreme)


# %%
def decadal_extrc_tsurf(index: xr.DataArray, temp: xr.DataArray, hlayers = None):
    """
    extract the extreme count and the mean surface temperature every ten years.
    **Arguments**
        *index* the DataArray of pc index
        *temp* the -fldmean, -yearmean -ensmean surface temperauter.
    **Return**
        *ext_counts* the DataArray of extreme count, *time, *extr_type, *mode
        *t_surf_mean* the mean t_surface,
        the time here use the first year of the decade.
    """
    if hlayers is not None:
        index = index.sel(hlayers=hlayers)

    # tsurf increase
    temp = temp - temp.isel(time = 0)

    ext_counts = []
    t_surf_mean = []
    for i in range(0, index.time.size, 10):

        period = slice(i, i + 10)

        period_pc = index.isel(time=period)
        period_tm = temp.isel(time=period)

        # extreme count
        period_ext_count = extreme.extreme_count_xr(period_pc)
        period_mean_t = period_tm.mean(dim="time")

        ext_counts.append(period_ext_count)
        t_surf_mean.append(period_mean_t)

    periods = index.isel(time=np.arange(0, index.time.size, 10)).time
    ext_counts = xr.concat(ext_counts, dim=periods)
    t_surf_mean = xr.concat(t_surf_mean, dim=periods)
    return ext_counts, t_surf_mean

#%%
def extCount_tsurf_scatter(ext_counts, t_surf, hlayers = None):
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
        xlim=(-1, 6),
        suptitle=f"extreme counts v.s surface temperature in decadal time scale",
        xlabel="temperature (K)",
        ylabel="extreme count",
        grid=False,
        leftlabels=["NAO", "EA"],
        toplabels=["pos", "neg"],
        xminorticks="null",
        yminorticks="null",
        ylim = (0,55),
    )
    if hlayers is not None:
        ext_counts = ext_counts.sel(hlayers=hlayers)

    for j, extr_type in enumerate(ext_counts.extr_type):
        for i, mode in enumerate(ext_counts.mode):

            # true values
            true = ext_counts.sel(extr_type=extr_type, mode=mode,confidence="true")
            low =  ext_counts.sel(extr_type=extr_type, mode=mode,confidence="low")
            high = ext_counts.sel(extr_type=extr_type, mode=mode,confidence="high")

            axes[i, j].errorbar(
                x=t_surf,
                y=true,
                yerr=[(true-low), (high-true)],
                fmt='o', linewidth=2, capsize=6)

# %%
