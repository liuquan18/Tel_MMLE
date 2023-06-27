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

# %%
def decadal_extrc_tsurf(index: xr.DataArray, ext_counts_xr = None, temp: xr.DataArray = None, plev=None,ci = 'AR1'):
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



    # start time
    time_s = index.time[::10]
    # end time
    time_e = index.time[9::10]

    # create slice for each decade
    decade_slice = [slice(s, e) for s, e in zip(time_s, time_e)]

    ext_counts = []
    t_surf_mean = []


    for time in decade_slice:
        print(f" extreme counting in the decade of {time.start.dt.year.values} - {time.stop.dt.year.values}")

        period_pc = index.sel(time=time)
        # ensure that there are 10 years of data in period_pc
        if period_pc.time.size != 10:
            print(f" the length of the period is {len(period_pc.time)}, skip this period")
            # rasing a warning
            warnings.warn(f" the length of the period is {len(period_pc.time)}")
            break
        time_tag = period_pc.time[0] # for reference 

        # extreme count
        if ext_counts_xr is None:
            period_ext_count = extreme.extreme_count_xr(period_pc, ci=ci)
            period_ext_count['time'] = time_tag
            # set time as the new dimension
            period_ext_count = period_ext_count.expand_dims('time')
            ext_counts.append(period_ext_count)

        # tsurf
        if temp is not None:
            period_tm = temp.sel(time = period_pc.time, method = 'nearest')
            period_mean_t = period_tm.mean(dim="time")
            period_mean_t['time'] = time_tag
            period_mean_t = period_mean_t.expand_dims('time')
            t_surf_mean.append(period_mean_t)

    if ext_counts_xr is None:
        ext_counts_xr = xr.concat(ext_counts, dim='time')
    if temp is not None:
        t_surf_mean_xr = xr.concat(t_surf_mean, dim='time')
        
    return ext_counts_xr, t_surf_mean_xr
