#%%
import xarray as xr
import numpy as np
import src.extreme.extreme_ci as extreme
import proplot as pplt
import warnings
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import statsmodels.api as sm

#%%

################### extreme event count v.s surface temperature scatter ###################
#%%
def extCount_tsurf_scatter(
    ext_counts, t_surf, ylim=(0, 55), xlim=(-1, 5), xlabel="temperature (K)"
):
    """
    plot the scatter plot of extreme event count v.s surface temperature
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

    _scatter_extrcVStsurf(ext_counts, t_surf, axes)


def _scatter_extrcVStsurf(ext_counts, t_surf, axes):
    t_surf = t_surf.sel(time=ext_counts.time, method="nearest")
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


#################### extreme event count profile ####################
#%%
# plot function, x-axis is the extreme events count, y-axis is the pressure level
# vertical line, and fill_betweenx the confidence interval
def _plot_extreme_count(ext_count, ax=None, label=None, colored=False):
    """
    plot the vertical profile of extreme event count for a single mode and extreme type
    """
    color = None
    style = None
    if colored:
        if label == "first10":
            color = "#1f77b4"
            style = color
        elif label == "last10":
            color = "#ff7f0e"
            style = color
    else:
        if label == "first10":
            color = "gray"
            style = "k-"
        elif label == "last10":
            color = "gray"
            style = "k--"

    if ax is None:
        ax = plt.gca()

    y = ext_count.plev / 100
    true = ext_count.sel(confidence="true").values
    low = ext_count.sel(confidence="low").values
    high = ext_count.sel(confidence="high").values

    # plot the confidence interval
    ax.fill_betweenx(y, low, high, color=color, alpha=0.3)
    # plot the true value
    line = ax.plot(true, y, style, linewidth=1, label=label)
    return line


def extreme_count_profile(first_count, last_count, colored=False, **kwargs):
    """
    plot the extreme event count profile for the NAO and EA,
    and positive and negative extreme events
    """
    # parameters from kwargs
    xlim = kwargs.pop("xlim", None)

    fig = pplt.figure(
        # space=0,
        refwidth="20em",
    )
    axes = fig.subplots(nrows=2, ncols=2)
    axes.format(
        abc=True,
        abcloc="ul",
        abcstyle="a",
        xlabel="Extreme event count",
        ylabel="Pressure level (hPa)",
        suptitle="Extreme event count profile",
        ylim=(1000, 200),
        xminorticks="null",
        yminorticks="null",
        grid=False,
        toplabels=("pos", "neg"),
        leftlabels=("NAO", "EA"),
        xlocator=20,
    )

    # plot the extreme event count profile
    labels = ["first10", "last10"]
    # the default color of matplotlib

    for i, extreme_count in enumerate([first_count, last_count]):
        _plot_extreme_count(
            extreme_count.sel(mode="NAO", extr_type="pos"),
            axes[0, 0],
            label=labels[i],
            colored=colored,
        )
        _plot_extreme_count(
            extreme_count.sel(mode="NAO", extr_type="neg"),
            axes[0, 1],
            label=labels[i],
            colored=colored,
        )

        _plot_extreme_count(
            extreme_count.sel(mode="EA", extr_type="pos"),
            axes[1, 0],
            label=labels[i],
            colored=colored,
        )
        _plot_extreme_count(
            extreme_count.sel(mode="EA", extr_type="neg"),
            axes[1, 1],
            label=labels[i],
            colored=colored,
        )
    for ax in axes:
        ax.set_xlim(xlim)
    # add legend
    axes[0, 0].legend(loc="lr", ncols=1, frame=True)


########################## MMLEA slope ################################
# %%
# calcualte slope
def calc_slope( extreme_count,tsurf):
    if tsurf is not None:
        x = tsurf.squeeze().values
    else:
        x = np.arange(len(extreme_count.time))
    y = extreme_count.sel(confidence="true").pc.values

    model = sm.OLS(y, sm.add_constant(x)).fit()
    slope = model.params[1]
    conf_int = model.conf_int()[1]
    return slope, conf_int


def slope_err(extr, tsurf):

    slope_NAO, conf_int_NAO = calc_slope(extr.sel(mode="NAO"), tsurf = None)
    slope_EA, conf_int_EA = calc_slope(extr.sel(mode="EA"), tsurf = None)

    yerr = np.array([[slope_NAO - conf_int_NAO[0]], [conf_int_NAO[1] - slope_NAO]])
    xerr = np.array([[slope_EA - conf_int_EA[0]], [conf_int_EA[1] - slope_EA]])
    return slope_NAO, slope_EA, yerr, xerr


def slope_models(
    extrs, tsurfs, models, axs, colors, ensemble_size=None, alpha=0.7, time="all"
):
    """
    plot the slope of extreme event count profile for all models
    """

    extr_types = ["pos", "neg"]
    modes = ["NAO", "EA"]
    if len(colors) != len(models):
        colors = [colors] * len(models)

    for i, model in enumerate(models):

        for j, extr_type in enumerate(extr_types):
            ax = axs[j]

            # select time period
            if time != "all" and model != "MPI_GE_onepct":
                time = np.datetime64(time)
                extr = extrs[model].sel(extr_type=extr_type, time=slice(time, None))
                tsurf = tsurfs[model].sel(time=extr.time, method="nearest")
            else:
                extr = extrs[model].sel(extr_type=extr_type)
                tsurf = tsurfs[model]
            tsurf = tsurf - tsurf[0]  # increase

            slope_NAO, slope_EA, yerr, xerr = slope_err(extr, tsurf)

            # calculate size of circle based on ensemble size
            if ensemble_size is None:
                size = 7
            else:
                size = ensemble_size[i] / 4

            scatter = ax.errorbar(
                slope_NAO,
                slope_EA,
                yerr=yerr,
                xerr=xerr,
                fmt="o",
                capsize=3,
                label=model,
                color=colors[i],
                markersize=size,
                alpha=alpha,
                markeredgewidth=0.5,
                # no edge color
                markeredgecolor="none",
            )
            ax.set_xlabel(f"Slope ({modes[0]})")
            ax.set_ylabel(f"Slope ({modes[1]})")
            ax.axhline(y=0, color="k", linewidth=0.5)
            ax.axvline(x=0, color="k", linewidth=0.5)
    return axs


def handle_label_models(colors, models):
    lines = [mlines.Line2D([], [], color=c, marker="o", markersize=5) for c in colors]
    labels = models
    return lines, labels


# Create a scatter plot of the slopes for each model and extreme type
def mmle_slope_scatter(
    extrs, tsurfs, extrs_rand, tsurfs_rand, tsurf="ens_fld_year_mean", time="all"
):
    """
    plot the slope of extreme event count profile for multiple models
    """

    fig, axs = pplt.subplots(nrows=2, ncols=2, sharex=False, sharey=False)

    axs.format(
        suptitle="Slopes of extreme counts vs. temperature",
        abc=True,
        grid=False,
        xtickminor=False,
        ytickminor=False,
    )
    models = [
        "MPI_GE_onepct",
        "MPI_GE",
        "CanESM2",
        "CESM1_CAM5",
        "MK36",
        "GFDL_CM3",
    ]
    colors_model = ["tab:red", "C1", "tab:blue", "tab:purple", "tab:cyan", "C4"]
    ensemble_size = [100, 100, 50, 40, 30, 20]

    # Get a list of nine evenly spaced colors from the colormap
    ens_size = np.arange(20, 101, 10)

    slope_models(
        extrs,
        tsurfs,
        models,
        axs[0, :],
        colors_model,
        ensemble_size=ensemble_size,
        time=time,
    )
    slope_models(
        extrs_rand,
        tsurfs_rand,
        np.arange(20, 101, 10),
        axs[1, :],
        "tab:grey",
        ensemble_size=ens_size,
        alpha=0.5,
    )
    # legend
    # Add the legend with the custom handler
    handles_color, labels_model = handle_label_models(colors_model, models)

    axs[0, 1].legend(
        handles_color,
        labels_model,
        loc="b",
        ncols=3,
        frame=False,
        facecolor="none",
        bbox_to_anchor=(-0.1, 0.6),
        # space from the plot
        columnspacing=5,
    )

    axs[1, 0].legend(
        loc="ll",
        ncols=1,
        facecolor="none",
        # move the legend higher
        bbox_to_anchor=(0, 0.25, 1, 1),
        frame=False,
        title="ens size",
        # the row space
    )

    axs[:, 0].format(
        xlim=(-5.1, 7.9),
    )

    axs[0, :].format(ylim=(-5.2, 10.7))

    axs[:, 1].format(
        xlim=(-4.2, 9.2),
    )

    axs[1, :].format(ylim=(-2.4, 10.4))

    if tsurf == "tropical_arctic_gradient":
        axs[:, 0].format(
            xlim=(-5.1, 2.1),
        )

        axs[0, :].format(ylim=(-5.2, 2.7))

        axs[:, 1].format(
            xlim=(-5.1, 2.1),
        )

        axs[1, :].format(ylim=(-5.2, 2.7))


####################### MMLE slope scatter with different tsurf #############################
#%%
def slope_diff_tsurf(extrs, tsurf_gmst, NA_tsurf, tropical_arctic_gradient, time="all"):
    """
    plot the slope of extreme event count profile for multiple models, and different tsurf as the x axis
    """
    fig, axs = pplt.subplots(nrows=3, ncols=2, sharex=False, sharey=True)
    axs.format(
        suptitle="Slopes of extreme counts vs. temperature",
        abc=True,
        grid=False,
        xtickminor=False,
        ytickminor=False,
        leftlabels=["GMST", "NA", "Tropical-Arctic"],
        toplabels=["pos", "neg"],
    )
    models = [
        "MPI_GE_onepct",
        "MPI_GE",
        "CanESM2",
        "CESM1_CAM5",
        "MK36",
        "GFDL_CM3",
    ]
    colors_model = ["tab:red", "C1", "tab:blue", "tab:purple", "tab:cyan", "C4"]
    ensemble_size = [100, 100, 50, 40, 30, 20]

    slope_models(
        extrs,
        tsurf_gmst,
        models,
        axs[0, :],
        colors_model,
        ensemble_size=ensemble_size,
        time=time,
    )
    slope_models(
        extrs,
        NA_tsurf,
        models,
        axs[1, :],
        colors_model,
        ensemble_size=ensemble_size,
        time=time,
    )
    slope_models(
        extrs,
        tropical_arctic_gradient,
        models,
        axs[2, :],
        colors_model,
        ensemble_size=ensemble_size,
        time=time,
    )

    handles_color, labels_model = handle_label_models(colors_model, models)

    fig.legend(
        handles_color,
        labels_model,
        loc="b",
        ncols=3,
        frame=False,
        facecolor="none",
        # bbox_to_anchor=(-0.1, 0.6),
        # space from the plot
        columnspacing=5,
    )

    axs[2, :].format(
        # reverse the x and y axis
        xlim=(1.9, -3),
        ylim=(1.9, -3.5),
    )


############################### MMLE line plot #########################################
#%%
def mmle_line_plot(
    extrs, tsurfs, extrs_rands, tsurfs_rands, tsurf="ens_fld_year_mean", time="all"
):
    params = {
        "ytick.color": "w",
        "xtick.color": "w",
        "axes.labelcolor": "w",
        "axes.edgecolor": "w",
        "tick.labelcolor": "w",
        "text.color": "w",
        "font.size": 20,
    }

    pplt.rc.update(params)

    fig, axs = pplt.subplots(
        nrows=1,
        ncols=2,
        sharex=True,
        sharey=False,
        figsize=(18, 12),
        facecolor="k",
        wspace=8,
        # dpi = 300,
    )

    axs.format(
        # suptitle="extreme event occurence vs. global mean temperature increase",
        abc=True,
        grid=False,
        xtickminor=False,
        ytickminor=False,
        xlabel="global mean temperature increase (K)",
        fontsize=25,
        toplabels=["pos", "neg"],
        facecolor="k",
        suptitle_kw = dict(color='w')
    )
    models = [
        "MPI_GE_onepct",
        "MPI_GE",
        "CanESM2",
        "CESM1_CAM5",
        "MK36",
        "GFDL_CM3",
    ]
    colors_model = ["tab:red", "C1", "tab:purple", "tab:blue", "tab:green", "C4"]
    # colors_model = ['r','#A7554C','#A776A3','#44ADC8','#5FD18D','#F4D653']
    ensemble_size = [100, 100, 50, 40, 30, 20]

    # rows for different modes, columns for pos and neg, colors for different models
    for r, mode in enumerate(["NAO"]):
        for c, extr_type in enumerate(["pos", "neg"]):
            for i, model in enumerate(models):

                # the ensemble size for each model
                ens_size = ensemble_size[i]

                # single data
                if time != "all" and model != "MPI_GE_onepct":
                    time = np.datetime64(time)
                    extrc = extrs[model].sel(
                        mode=mode, extr_type=extr_type, time=slice(time, '2091-01-01')
                    )
                    tsurf = tsurfs[model].sel(time=extrc.time, method="nearest")
                else:
                    extrc = extrs[model].sel(mode=mode, extr_type=extr_type)
                    tsurf = tsurfs[model].sel(time=extrc.time, method="nearest")

                extrc_rand = extrs_rands[ens_size].sel(mode=mode, extr_type=extr_type)
                tsurf_rand = tsurfs_rands[ens_size].sel(
                    time=extrc_rand.time, method="nearest"
                )

                # plot the line
                im = line_single(
                    extrc,
                    tsurf,
                    ensemble_size[i],
                    axs[r, c],
                    color=colors_model[i],
                    label=f"{model} ({str(ens_size)})",
                )
    axs[0].format(
        ylabel="extreme occurence",
        ylim=(0, 118),
        ylocator=20,
    )
    axs[1].format(
        ylim=(0, 98),
        # yticks every 20
        ylocator=20,
    )
    axs[1].legend(
        loc="r",
        ncols=1,
        fontsize=40,
        columnspacing=2,
        labelspacing=4,
        facecolor='none',
        frame=False,
        borderaxespad=5,
        bbox_to_anchor=(1.1, 0.5),
    )

    for ax in axs:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    return fig


def line_single(
    extrc, tsurf, ens_size, ax, extrc_rand=None, tsurf_rand=None, color="k", label=None
):
    """
    line for just one dataset.
    """
    tsurf = tsurf.values
    tsurf = tsurf - tsurf[0]
    extrc = extrc.sel(confidence="true").pc.squeeze().values

    ax.line(
        x=tsurf,
        y=extrc,
        marker="o",
        color=color,
        markersize=6,
        label=label,
        linewidth=4,
    )

    if extrc_rand is not None:
        tsurf_rand = tsurf_rand.values
        tsurf_rand = tsurf_rand - tsurf_rand[0]
        extr_rand = extrc_rand.sel(confidence="true").pc.squeeze().values

        ax.plot(
            x=tsurf_rand,
            y=extr_rand,
            marker="o",
            color=color,
            alpha=0.5,
            markersize=6,
            linestyle="--",
        )
    else:
        pass
    ax.set_xlim(-1, 5.6)


#%%
def mmle_line_sum_plot(
    extrs, tsurfs, extrs_rands, tsurfs_rands, tsurf="ens_fld_year_mean", time="all"
):

    fig, axs = pplt.subplots(
        nrows=1, ncols=2, sharex=False, sharey=False, figsize=(10, 20)
    )

    axs.format(
        suptitle="Slopes of extreme counts vs. temperature",
        abc=True,
        grid=False,
        xtickminor=False,
        ytickminor=False,
    )
    models = [
        "MPI_GE_onepct",
        "MPI_GE",
        "CanESM2",
        "CESM1_CAM5",
        "MK36",
        "GFDL_CM3",
    ]
    # colors_model = ["tab:red", "C1", "tab:blue", "tab:purple", "tab:cyan", "C4"]
    colors_model = ["tab:black", "#A7554C", "#A776A3", "#44ADC8", "#5FD18D", "#F4D653"]
    ensemble_size = [100, 100, 50, 40, 30, 20]

    # rows for different modes, columns for pos and neg, colors for different models
    for c, extr_type in enumerate(["pos", "neg"]):
        for i, model in enumerate(models):

            # the ensemble size for each model
            ens_size = ensemble_size[i]

            # single data
            if time != "all" and model != "MPI_GE_onepct":
                time = np.datetime64(time)
                extrc = (
                    extrs[model]
                    .sum(dim="mode")
                    .sel(extr_type=extr_type, time=slice(time, None))
                )
                tsurf = tsurfs[model].sel(time=extrc.time, method="nearest")
            else:
                extrc = extrs[model].sum(dim="mode").sel(extr_type=extr_type)
                tsurf = tsurfs[model].sel(time=extrc.time, method="nearest")

            extrc_rand = extrs_rands[ens_size].sum(dim="mode").sel(extr_type=extr_type)
            tsurf_rand = tsurfs_rands[ens_size].sel(
                time=extrc_rand.time, method="nearest"
            )

            # plot the line
            line_single(
                extrc,
                tsurf,
                ensemble_size[i],
                axs[c],
                # extrc_rand=extrc_rand,
                # tsurf_rand=tsurf_rand,
                color=colors_model[i],
            )
