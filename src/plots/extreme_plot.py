#%%
import xarray as xr
import numpy as np
import src.extreme.extreme_ci as extreme
import proplot as pplt
import warnings
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import statsmodels.api as sm

import matplotlib.patches as patches
import pandas as pd

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
    params = {
        "ytick.color": "k",
        "xtick.color": "k",
        "axes.labelcolor": "k",
        "axes.edgecolor": "k",
        "tick.labelcolor": "k",
        "text.color": "k",
        "font.size": 20,
    }

    pplt.rc.update(params)

    fig, axes = pplt.subplots(nrows=1, ncols=2, figwidth=8, span=False, sharey = True, sharex = True,facecolor="black")

    axes.format(
        abc="a",
        abcloc="ul",
        # xlim=xlim,
        suptitle=f"extreme counts v.s surface temperature in decadal time scale",
        xlabel=xlabel,
        ylabel="extreme count",
        grid=False,
        # leftlabels=["NAO", "EA"],
        toplabels=["pos", "neg"],
        xminorticks="null",
        yminorticks="null",
        facecolor="white",
        # ylim=ylim,
    )

    t_surf = t_surf.sel(time=ext_counts.time, method="nearest")
    t_surf = t_surf - t_surf[0]
    for i, mode in enumerate(['NAO']):
        for j, extr_type in enumerate(ext_counts.extr_type):
            # data preparation
            true = ext_counts.sel(extr_type=extr_type, mode=mode, confidence="true")
            low = ext_counts.sel(extr_type=extr_type, mode=mode, confidence="low")
            high = ext_counts.sel(extr_type=extr_type, mode=mode, confidence="high")

            # for the data with plev

            t = t_surf

            true = true[[0,-1]]
            t = t[[0,-1]]
            low = low[[0,-1]]
            high = high[[0,-1]]

            axes[i, j].errorbar(
                x=t,
                y=true,
                yerr=[(true - low), (high - true)],
                fmt="o",
                linewidth=2,
                capsize=6,
                color = 'red',
            )

            axes[i, j].set_xlim(-1, 5)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

#################### extreme event count profile ####################
#%%
# plot function, x-axis is the extreme events count, y-axis is the pressure level
# vertical line, and fill_betweenx the confidence interval
def plot_extreme_count(ext_count, ax=None, label=None, conf_color = 'gray',line_color = 'k', linewidth = 2):
    """
    plot the vertical profile of extreme event count for a single mode and extreme type
    """
    color = None
    style = None
    if 'first' in label:
        color = conf_color
        style = f"{line_color}-"
    elif 'last' in label:
        color = conf_color
        style = f"{line_color}--"

    if ax is None:
        ax = plt.gca()

    y = ext_count.plev / 100
    if y.min() < 200:
        y = ext_count.plev
    true = ext_count.sel(confidence="true").values
    low = ext_count.sel(confidence="low").values
    high = ext_count.sel(confidence="high").values

    # plot the confidence interval
    ax.fill_betweenx(y, low, high, color=color, alpha=0.3)
    # plot the true value
    line = ax.plot(true, y, style, linewidth=linewidth, label=label)
    return line


def extreme_count_profile(first_count, last_count, colored=False, **kwargs):
    """
    plot the extreme event count profile for the NAO and EA,
    and positive and negative extreme events
    """
    # parameters from kwargs
    xlim = kwargs.pop("xlim", None)
    axes = kwargs.pop("axes", None)
    if axes is None:
        # build the figure
        fig = pplt.figure(
            # space=0,
            refwidth="20em",
            facecolor="white",
        )
        axes = fig.subplots(nrows=1, ncols=2)
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
            # leftlabels=("NAO", "EA"),
            # xlocator=20,
        )

    # plot the extreme event count profile
    labels = ["first10", "last10"]
    # the default color of matplotlib

    for i, extreme_count in enumerate([first_count, last_count]):
        plot_extreme_count(
            extreme_count.sel(mode="NAO", extr_type="pos"),
            axes[0, 0],
            label=labels[i],
            colored=colored,
        )
        plot_extreme_count(
            extreme_count.sel(mode="NAO", extr_type="neg"),
            axes[0, 1],
            label=labels[i],
            colored=colored,
        )

        # _plot_extreme_count(
        #     extreme_count.sel(mode="EA", extr_type="pos"),
        #     axes[1, 0],
        #     label=labels[i],
        #     colored=colored,
        # )
        # _plot_extreme_count(
        #     extreme_count.sel(mode="EA", extr_type="neg"),
        #     axes[1, 1],
        #     label=labels[i],
        #     colored=colored,
        # )
    for ax in axes:
        ax.set_xlim(xlim)
    # add legend
    axes[0, 0].legend(loc="lr", ncols=1, frame=True)

    return axes


########################## MMLEA slope ################################
# %%
# calcualte slope
def calc_slope( extreme_count,tsurf):
    if tsurf is not None:
        x = tsurf.squeeze().values
    else:
        x = np.arange(len(extreme_count.time))
    try:
        y = extreme_count.sel(confidence="true").values
    except KeyError:
        y = extreme_count.values

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
    extrs, tsurfs, extrs_rand, tsurfs_rand, tsurf="ens_fld_year_mean", time="all",
    models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"],
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
    models_all = [
        "MPI_GE_onepct",
        "MPI_GE",
        "CanESM2",
        "CESM1_CAM5",
        "MK36",
        "GFDL_CM3",
    ]
    colors_model = ["tab:red", "C1", "tab:blue", "tab:purple", "tab:cyan", "C4"]
    model_color = dict(zip(models_all, colors_model))
    ensemble_size = [100, 100, 50, 40, 30, 20]

    # Get a list of nine evenly spaced colors from the colormap
    ens_size = np.arange(20, 101, 10)
    extrs = {k: extrs[k] for k in models}
    tsurfs = {k: tsurfs[k] for k in models} # select the subset of models
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
def mmle_tsurf_line(
    extrs, tsurfs, extrs_rands, tsurfs_rands, tsurf="ens_fld_year_mean", time="all",x_var = 'time'
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
        sharey=True,
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
                elif model == "MPI_GE_onepct":
                    extrc = extrs[model].sel(mode=mode, extr_type=extr_type)
                    tsurf = tsurfs[model].sel(time=extrc.time, method="nearest")
                    tsurf['time'] = pd.date_range(time, freq='10Y',periods = 14)
                else:
                    extrc = extrs[model].sel(mode=mode, extr_type=extr_type)
                    tsurf = tsurfs[model].sel(time=extrc.time, method="nearest")

                extrc_rand = extrs_rands[ens_size].sel(mode=mode, extr_type=extr_type)
                tsurf_rand = tsurfs_rands[ens_size].sel(
                    time=extrc_rand.time, method="nearest"
                )

                # plot the line
                im = extrc_tsurf_line_single(
                    extrc,
                    tsurf,
                    axs[r, c],
                    color=colors_model[i],
                    label=f"{model} ({str(ens_size)})",
                    x_var=x_var,
                )
    axs[0].format(
        ylabel="extreme occurence",
        ylim=(0, 110),
        ylocator=20,
        xlocator=1,
    )
    axs[1].format(
        ylim=(0, 110),
        # yticks every 20
        ylocator=20,
        # xlocator=1,
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
        ax.set_yticks([0, 20, 40, 60, 80, 100])
    return fig


def extrc_tsurf_line_single(
    extrc, tsurf, ax, extrc_rand=None, tsurf_rand=None, color="k", label=None,x_var = 'time'
):
    """
    line for just one dataset.
    """
    tsurf = tsurf
    tsurf = tsurf - tsurf[0]
    extrc = extrc.sel(confidence="true").pc.squeeze().values

    if x_var == 'tsurf':
        x = tsurf.values
    elif x_var == 'time':
        x = tsurf.time.dt.year.values
    ax.line(
        x = x,
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

def extrc_time_line(extrcs, **kwargs):
    ylim = kwargs.pop("ylim", (20, 280))

    gs = pplt.GridSpec(nrows=1, ncols=2)
    fig = pplt.figure(refwidth=2.2, refheight = 5.2, span=False, share="labels")
    # the right order of the models
    models_legend = [
    "MPI_GE_onepct (100)",
    "MPI-GE (100)",
    "CanESM2 (50)",
    "CESM1-CAM5 (40)",
    "MK3.6 (30)",
    "GFDL-CM3 (20)",
]

    lines = []
    for r, mode in enumerate(['NAO']):
        for c, extr_type in enumerate(['pos','neg']):
            ax = fig.subplot(gs[r, c])
            line = extrc_time_line_single(extrcs,  extr_type, ax,ylim = ylim)
            lines.append(line)
    fig.legend(
    lines,
    labels=models_legend,
    ncols=3,
    loc="b",
)

    return fig

def extrc_time_line_single(extrcs, extr_type, ax, ylim = (20, 280),mode = 'NAO',ci = False,
                               models = ["MPI_GE_onepct","MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
):
    models_all = ["MPI_GE_onepct","MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3", "MPI_GE_RCP45"]
    colors_model = ["red", "C1", "tab:purple", "tab:blue", "tab:green", "C4", "k"]
    model_color = dict(zip(models_all, colors_model))

    lines = []
    for model in models:
        try:
            extrc = extrcs[models.index(model)]
        except KeyError:
            extrc = extrcs[model]
        line = extrc.sel(extr_type=extr_type,mode = mode,confidence = 'true').plot.line(
                ax=ax, 
                label=model_color[model],
                x = 'time',
                color = model_color[model],
                linewidth = 2,)
        
        if ci:
            # fill between the confidence interval ['low','high']
            ax.fill_between(
                extrc.time,
                extrc.sel(extr_type=extr_type,mode = mode,confidence = 'low').values,
                extrc.sel(extr_type=extr_type,mode = mode,confidence = 'high').values,
                color = model_color[model],
                alpha = 0.15,
            )
        lines.append(line)

            
    ax.format(
            ylim=ylim,
            ylabel="Extreme counts",
            xlabel="Year",
            suptitle="",
            titleloc="uc",
            ylocator=20,
            yminorlocator="null",
            grid=False,
            title = '',
        )
    return ax, lines



#%%
# for reananlysis data

def reananlysis_bar(first_extrc, first_err, last_extrc, last_err, ax, 
                    x = [0.2,0.8],width = 0.4,facecolor = 'none',edgecolor = 'black',linewidth = 1,errcolor = 'black'):

    ax.bar(
        x =x,
        height = [first_extrc,last_extrc],
        width = width,
        edgecolor = edgecolor,
        facecolor = facecolor,
        linewidth = 1,
        align = 'center',
        zorder = 9,
    )
    ax.errorbar(
        x = x,
        y = [first_extrc,last_extrc],
        yerr = [first_err,last_err],
        color = errcolor,
        linewidth = 2,
        fmt='none',
        zorder = 10,
    )
    return ax
#%%

#%%    
def format_ens_size(ens_size,tick_number):
    if ens_size == 70:
        formater = '100'
    else:
        formater = str(ens_size )
    return formater

def plot_errorbar(x,slope,low,high,color,ax,width = 2):
    line = ax.hlines(x1 = x-width,x2 = x,y = slope,color = color,linewidth = 2)
    x1 = x-width
    x2 = x
    y1 = low
    y2 = high
    rect = patches.Rectangle((x1, y1), x2-x1, y2-y1, linewidth=1, edgecolor='none', facecolor=color,alpha = 0.5)
    bar = ax.add_patch(rect)
    return line,bar

def plot_unfill_errbar(x,slope,low,high,color,ax,width = 2):
    line = ax.hlines(x1 = x-width,x2 = x,y = slope,color = color,linewidth = 2)
    x1 = x-width
    x2 = x
    y1 = low
    y2 = high
    rect_out = patches.Rectangle((x1, y1), x2-x1, y2-y1, linewidth=1, edgecolor=color, facecolor='none',alpha = 0.5)
    bar = ax.add_patch(rect_out)
    return line,bar


#%%
def extrc_slope_line(slopes,ax,mode = 'NAO',extr_type = 'pos',
                     models = ["MPI_GE_onepct", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"],rand = False,shift = 0):
    
    models_all = ["MPI_GE_onepct","MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    colors_model = ["red", "C1", "tab:purple", "tab:blue", "tab:green", "C4"]
    ens_size = [70, 70, 50, 40, 30, 20]  # 70 for the position of 100
    model_color = dict(zip(models_all, colors_model))
    model_size = dict(zip(models_all, ens_size))

    lines = []
    bars = []
    if not rand:
        for model in models[::-1]: # from low to high
            # select time period
            slope = slopes[model].sel(extr_type = extr_type,mode = mode,slopes = 'true').values
            low = slopes[model].sel(extr_type = extr_type,mode = mode,slopes = 'low').values
            high = slopes[model].sel(extr_type = extr_type,mode = mode,slopes = 'high').values

            # a shift for the x position
            # define the shift as 0 when rand = True, and 1 when rand = False
            model_size['MPI_GE'] = 70+shift
            model_size['MPI_GE_onepct'] = 70-shift

            line,bar = plot_errorbar(
                x = model_size[model]+shift,
                slope = slope,
                low = low,
                high = high,
                color = model_color[model],
                ax = ax,
            )
            lines.append(line)

            bars.append(bar)

    else:
        for ens_size in [20,30,40,50,100]:
            # select time period
            slope = slopes[ens_size].sel(extr_type = extr_type,mode = mode,slopes = 'true').values
            low = slopes[ens_size].sel(extr_type = extr_type,mode = mode,slopes = 'low').values
            high = slopes[ens_size].sel(extr_type = extr_type,mode = mode,slopes = 'high').values

            # the position of 100 is 70
            if ens_size == 100:
                x = 70 + 2
                line,bar = plot_errorbar(
                    x = x,
                    slope = slope,
                    low = low,
                    high = high,
                    color = model_color['MPI_GE'],
                    ax = ax,
                )
            else:
                x = ens_size + 2
                line,bar = plot_unfill_errbar(
                    x = x,
                    slope = slope,
                    low = low,
                    high = high,
                    color = model_color['MPI_GE'],
                    ax = ax,
                )
            


    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_ens_size))

    return ax,bars
# %%
