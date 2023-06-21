# %%
import proplot as pplt
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

# %%
def read_extreme_counts(plev = 50000,standard = 'first',tsurf = 'ens_fld_year_mean',fixed_pattern = 'decade'):
    models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
    # load data
    # different models
    extrs = {}
    tsurfs = {}
    for model in models:
        odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
        prefix = f"plev_{plev}_{fixed_pattern}_{standard}_"
        extrs[model] = xr.open_dataset(
            f"{odir}{prefix}extre_counts.nc"
        ).squeeze()
        tsurfs[model] = xr.open_dataset(
            f"{odir}{tsurf}.nc"
        )
        try:
            tsurfs[model] = tsurfs[model].tsurf
        except AttributeError:
            try:
                tsurfs[model] = tsurfs[model].ts
            except AttributeError:
                tsurfs[model] = tsurfs[model].tas

    # random sampled models
    extrs_rand = {}
    tsurfs_rand = {}
    for ens_size in np.arange(20, 101, 10):
        odir = f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/"
        extrs_rand[ens_size] = xr.open_dataset(
            f"{odir}plev_{plev}_{standard}_extreme_counts_{ens_size}.nc"
        ).squeeze()
        tsurfs_rand[ens_size] = xr.open_dataset(
            f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/{tsurf}.nc"
        ).tsurf

    return extrs, tsurfs, extrs_rand, tsurfs_rand 

# %%
# calcualte slope
def calc_slope(tsurf,extreme_counts, extr_type, mode):
    x = tsurf.values
    y = extreme_counts.sel(extr_type=extr_type, mode=mode, confidence="true").pc.values

    model = sm.OLS(y, sm.add_constant(x)).fit()
    slope = model.params[1]
    conf_int = model.conf_int()[1]
    return slope, conf_int


# %%
def plot_scatter(extrs,tsurfs, models, axs, colors, ensemble_size=None,alpha = 0.7,time = 'all'):
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
def plot_slope(plev = 50000,standard = 'first',tsurf = 'ens_fld_year_mean',time = 'all'):
    models = ["MPI_GE_onepct","MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
    extrs,tsurfs,extrs_rand,tsurfs_rand = read_extreme_counts(plev = plev,standard = standard,tsurf = tsurf)
    fig, axs = pplt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)

    axs.format(
        suptitle="Slopes of extreme counts vs. temperature",
        abc=True,
        grid=False,
        xtickminor=False,
        ytickminor=False,
    )

    colors_model = ["tab:red", "C1", "tab:blue", "tab:purple", "C4", "tab:cyan"]
    ensemble_size = [100, 100, 45, 40, 30, 20]

    # Get the "autumn" colormap
    cmap = plt.get_cmap("autumn")
    # Get a list of nine evenly spaced colors from the colormap
    ens_size = np.arange(20, 101, 10)

    plot_scatter(extrs,tsurfs, models, axs[0, :], colors_model, ensemble_size=ensemble_size,time = time)
    plot_scatter(extrs_rand,tsurfs_rand, np.arange(20, 101, 10), axs[1, :], "tab:red", ensemble_size=ens_size,alpha=0.5)
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

    plt.savefig(f'/work/mh0033/m300883/Tel_MMLE/docs/source/plots/MMLE/plev_{plev}_extrc_{tsurf}_{time}_slope.png')

#%%
def slope_diff_tsurf(plev = 50000,standard = 'first',tsurf = 'ens_fld_year_mean',time = 'all'):
    extrs,tsurf_gmst,_,_ = read_extreme_counts(plev = plev,standard = standard,tsurf = 'ens_fld_year_mean')
    extrs,NA_tsurf,_,_ = read_extreme_counts(plev = plev,standard = standard,tsurf = 'NA_tsurf')
    extrs,tropical_arctic_gradient,_,_ = read_extreme_counts(plev = plev,standard = standard,tsurf = 'tropical_arctic_gradient')

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
    tsurfs = ['ens_fld_year_mean','NA_tsurf','tropical_arctic_gradient']

    colors_model = ["tab:red", "C1", "tab:blue", "tab:purple", "C4", "tab:cyan"]
    ensemble_size = [100, 100, 45, 40, 30, 20]

    plot_scatter(extrs,tsurf_gmst, models, axs[0, :], colors_model, ensemble_size=ensemble_size,time = time)
    plot_scatter(extrs,NA_tsurf, models, axs[1, :], colors_model, ensemble_size=ensemble_size,time = time)
    plot_scatter(extrs,tropical_arctic_gradient, models, axs[2, :], colors_model, ensemble_size=ensemble_size,time = time)

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
    plt.savefig(f'/work/mh0033/m300883/Tel_MMLE/docs/source/plots/MMLE/plev_{plev}_extrc_tsurf_all_{time}_slope.png')

#%%




# %%
plot_slope(plev=50000,tsurf='ens_fld_year_mean',standard='first')
#%%
plot_slope(plev=50000,tsurf='NA_tsurf',standard='first')
#%%
plot_slope(plev=50000,tsurf='tropical_arctic_gradient',standard='first')

#%%
plot_slope(plev=50000,tsurf='ens_fld_year_mean',standard='first',time = '1950-01-01')
#%%
plot_slope(plev=50000,tsurf='NA_tsurf',standard='first',time = '1950-01-01')
#%%
plot_slope(plev=50000,tsurf='tropical_arctic_gradient',standard='first',time = '1950-01-01')
#%%

#%%
plot_slope(plev=50000,standard='temporal_ens')

# %%
plot_slope(plev=30000)

# %%
