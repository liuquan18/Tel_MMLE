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

    _scatter_extrcVStsurf(ext_counts, t_surf,  axes)


def _scatter_extrcVStsurf(ext_counts, t_surf,  axes):
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


# %%
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



# %%
# calcualte slope
def calc_slope(tsurf,extreme_counts, extr_type, mode):
    x = tsurf.squeeze().values
    y = extreme_counts.sel(extr_type=extr_type, mode=mode, confidence="true").pc.values

    model = sm.OLS(y, sm.add_constant(x)).fit()
    slope = model.params[1]
    conf_int = model.conf_int()[1]
    return slope, conf_int


# %%
def _slope_single(extrs,tsurfs, models, axs, colors, ensemble_size=None,alpha = 0.7,time = 'all'):
    """
    plot the slope of extreme event count profile for a single model
    """
    extr_types = ["pos", "neg"]
    modes = ["NAO", "EA"]
    if len(colors) != len(models):
        colors = [colors] * len(models)

    for i, extr_type in enumerate(extr_types):
        ax = axs[i]
        for j, model in enumerate(models):

            if time != 'all' and model != 'MPI_GE_onepct':
                time = np.datetime64(time)
                tsurfs[model] = tsurfs[model].sel(time = slice(time,None))
                extrs[model] = extrs[model].sel(time = slice(time,None))

            slope_NAO, conf_int_NAO = calc_slope(tsurfs[model],extrs[model], extr_type, "NAO")
            slope_EA, conf_int_EA = calc_slope(tsurfs[model],extrs[model], extr_type, "EA")

            yerr = np.array(
                [[slope_NAO - conf_int_NAO[0]], [conf_int_NAO[1] - slope_NAO]]
            )
            xerr = np.array([[slope_EA - conf_int_EA[0]], [conf_int_EA[1] - slope_EA]])

            # calculate size of circle based on ensemble size
            if ensemble_size is None:
                size = 7
            else:
                size = ensemble_size[j] / 4

            scatter = ax.errorbar(
                slope_NAO,
                slope_EA,
                yerr=yerr,
                xerr=xerr,
                fmt="o",
                capsize=3,
                label=model,
                color=colors[j],
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
# %%
def handle_label_models(colors,models):
    lines = [mlines.Line2D([], [], color=c, marker="o", markersize=5)for c in colors]
    labels = models
    return lines, labels


#%%
# Create a scatter plot of the slopes for each model and extreme type
def mmle_slope_scatter(extrs,tsurfs,extrs_rand,tsurfs_rand,tsurf = 'ens_fld_year_mean',time = 'all'):
    """
    plot the slope of extreme event count profile for multiple models
    """
    fig, axs = pplt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

    axs.format(
        suptitle="Slopes of extreme counts vs. temperature",
        abc=True,
        grid=False,
        xtickminor=False,
        ytickminor=False,
    )
    models = ["MPI_GE_onepct","MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
    colors_model = ["tab:red", "C1", "tab:blue", "tab:purple", "C4", "tab:cyan"]
    ensemble_size = [100, 100, 45, 40, 30, 20]

    # Get the "autumn" colormap
    cmap = plt.get_cmap("autumn")
    # Get a list of nine evenly spaced colors from the colormap
    ens_size = np.arange(20, 101, 10)

    _slope_single(extrs,tsurfs, models, axs[0, :], colors_model, ensemble_size=ensemble_size,time = time)
    _slope_single(extrs_rand,tsurfs_rand, np.arange(20, 101, 10), axs[1, :], "tab:grey", ensemble_size=ens_size,alpha=0.5)
    # legend 
    # Add the legend with the custom handler
    handles_color, labels_model = handle_label_models(colors_model,models)

    axs[0,1].legend(
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

    axs[1,0].legend(
        loc = "ll",
        ncols = 1,
        facecolor="none",
        # move the legend higher
        bbox_to_anchor=(0, 0.25, 1, 1),
        frame=False,
        title = 'ens size',
        # the row space
    )

    axs[:,0].format(
        xlim = (-5.1,7.9),
    )

    axs[0,:].format(
        ylim = (-5.2,10.7)
    )

    axs[:,1].format(
        xlim = (-4.2,9.2),
    )

    axs[1,:].format(
        ylim = (-2.4,10.4)
    )

    if tsurf == 'tropical_arctic_gradient':
        axs[:,0].format(
            xlim = (-5.1,2.1),
        )

        axs[0,:].format(
            ylim = (-5.2,2.7)
        )

        axs[:,1].format(
            xlim = (-5.1,2.1),
        )

        axs[1,:].format(
            ylim = (-5.2,2.7)
        )


#%%
def slope_diff_tsurf(extrs,tsurf_gmst,NA_tsurf,tropical_arctic_gradient,time = 'all'):
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
        leftlabels = ['GMST','NA','Tropical-Arctic'],
        toplabels = ['pos','neg'],
    )

    models = ["MPI_GE_onepct","MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
    colors_model = ["tab:red", "C1", "tab:blue", "tab:purple", "C4", "tab:cyan"]
    ensemble_size = [100, 100, 45, 40, 30, 20]

    _slope_single(extrs,tsurf_gmst, models, axs[0, :], colors_model, ensemble_size=ensemble_size,time = time)
    _slope_single(extrs,NA_tsurf, models, axs[1, :], colors_model, ensemble_size=ensemble_size,time = time)
    _slope_single(extrs,tropical_arctic_gradient, models, axs[2, :], colors_model, ensemble_size=ensemble_size,time = time)

    handles_color, labels_model = handle_label_models(colors_model,models)

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

    axs[2,:].format(
        # reverse the x and y axis
        xlim = (1.9,-3),
        ylim = (1.9,-3.5)
    )


