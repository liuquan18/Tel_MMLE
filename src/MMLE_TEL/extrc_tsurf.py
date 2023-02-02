import xarray as xr
import numpy as np
import src.extreme.period_pattern_extreme as extreme
import proplot as pplt

# %%
def decadal_extrc_tsurf(index: xr.DataArray, temp: xr.DataArray, hlayers: int = 50000):
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
        period_ext_count = extreme.period_extreme_count(period_pc)
        period_mean_t = period_tm.mean(dim="time")

        ext_counts.append(period_ext_count)
        t_surf_mean.append(period_mean_t)

    periods = index.isel(time=np.arange(0, index.time.size, 10)).time
    ext_counts = xr.concat(ext_counts, dim=periods)
    t_surf_mean = xr.concat(t_surf_mean, dim=periods)
    return ext_counts, t_surf_mean


def extCount_tsurf_scatter(ext_counts, t_surf):
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
        toplabels=["NAO", "EA"],
        leftlabels=["pos", "neg"],
    )

    extr_types = ["pos", "neg"]  # rows
    modes = ["NAO", "EA"]  # cols

    for i, extr_type in enumerate(ext_counts.extr_type):
        for j, mode in enumerate(ext_counts.mode):
            axes[i, j].scatter(
                x=t_surf,
                y=ext_counts.sel(extr_type=extr_type, mode=mode),
            )
