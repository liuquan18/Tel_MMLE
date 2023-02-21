#%%
import xarray as xr
import numpy as np
import src.extreme.extreme_ci as extreme
import proplot as pplt

# reimport extreme
import importlib

importlib.reload(extreme)


# %%
def decadal_extrc_tsurf(index: xr.DataArray, temp: xr.DataArray, plev=None):
    """
    extract the extreme count and the mean surface temperature every ten years.
    **Arguments**
        *index* the DataArray of pc index
        *temp* the -fldmean, -yearmean -ensmean surface temperauter.
    **Return**
        *ext_counts* the DataArray of extreme count, *time, *extr_type, *mode
        *t_surf_mean* the mean t_surface (the increase of the temperature)
        the time here use the first year of the decade.
    """
    if plev is not None:
        index = index.sel(plev=plev)

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


def corr_coef(ext_count, tsurf_increase):
    """calculate the correlation coefficient between extreme count and surface temperature increase"""
    func = lambda x: np.corrcoef(x, tsurf_increase.values)[0][1]
    return xr.apply_ufunc(
        func,
        ext_count,
        input_core_dims=[["time"]],
        output_core_dims=[[]],
        vectorize=True,
    )


#%%
def extCount_tsurf_scatter(ext_counts, t_surf, plev=None, ylim=(0, 55),xlim=(-1, 6),xlabel = 'temperature (K)'):
    """
    rows: pos/neg
    cols: NAO/EA
    scatter: extreme_count v.s surface temperature
    hue: different  dataset
    """
    r = corr_coef(ext_counts.sel(confidence="true"), t_surf)
    fig, axes = pplt.subplots(nrows=2, ncols=2, figwidth=8, span=False, share=False)

    axes.format(
        abc="a",
        abcloc="ul",
        xlim=xlim,
        suptitle=f"extreme counts v.s surface temperature in decadal time scale",
        xlabel=xlabel,
        ylabel="extreme count",
        grid=False,
        leftlabels=["NAO", "EA"],
        toplabels=["pos", "neg"],
        xminorticks="null",
        yminorticks="null",
        ylim=ylim,
    )
    if plev is not None:
        ext_counts = ext_counts.sel(plev=plev)
        r = r.sel(plev=plev)
        scatter_single_plev(ext_counts, t_surf, r, axes)

    else:
        for plev in ext_counts.plev:
            scatter_single_plev(ext_counts.sel(plev=plev), t_surf, r.sel(plev=plev), axes)

def scatter_single_plev(ext_counts, t_surf, r, axes):
    for j, extr_type in enumerate(ext_counts.extr_type):
        for i, mode in enumerate(ext_counts.mode):
            # true values
            true = ext_counts.sel(extr_type=extr_type, mode=mode, confidence="true")
            low = ext_counts.sel(extr_type=extr_type, mode=mode, confidence="low")
            high = ext_counts.sel(extr_type=extr_type, mode=mode, confidence="high")

            axes[i, j].errorbar(
                x=t_surf,
                y=true,
                yerr=[(true - low), (high - true)],
                fmt="o",
                linewidth=2,
                capsize=6,
            )

            axes[i, j].text(
                0.8,
                0.1,
                f"r={r.sel(extr_type=extr_type, mode=mode).values:.2f}",
                transform=axes[i, j].transAxes,
                fontsize="large",
            )


# %%
