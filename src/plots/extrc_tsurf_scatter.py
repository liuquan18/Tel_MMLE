#%%
import xarray as xr
import numpy as np
import src.extreme.extreme_ci as extreme
import proplot as pplt

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
    for j, extr_type in enumerate(ext_counts.extr_type):
        for i, mode in enumerate(ext_counts.mode):
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

# %%
def decadal_extrc_tsurf(index: xr.DataArray, temp: xr.DataArray, plev=None,ci = 'AR1'):
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
    print("counting the occureance of extremes and signifcance interval ...")
    if plev is not None:
        index = index.sel(plev=plev)

    # change the temp time into datetime64
    try:
        temp["time"] = temp.indexes["time"].to_datetimeindex()
    except AttributeError:
        pass

    ext_counts = []
    t_surf_mean = []
    for i in range(0, index.time.size, 10): # avoid the last period which is not 10 years
        
        period = slice(i, i + 10)

        period_pc = index.isel(time=period)
        period_tm = temp.sel(time = period_pc.time, method = 'nearest')
        time_tag = period_pc.time[0] # for reference 
        print(f" time: {time_tag.dt.year.values}")

        # extreme count
        period_ext_count = extreme.extreme_count_xr(period_pc, ci=ci)
        period_mean_t = period_tm.mean(dim="time")

        period_ext_count['time'] = time_tag
        period_mean_t['time'] = time_tag

        # set time as the new dimension
        period_ext_count = period_ext_count.expand_dims('time')
        period_mean_t = period_mean_t.expand_dims('time')

        ext_counts.append(period_ext_count)
        t_surf_mean.append(period_mean_t)

    ext_counts = xr.concat(ext_counts, dim='time')
    t_surf_mean = xr.concat(t_surf_mean, dim='time')
    return ext_counts, t_surf_mean
