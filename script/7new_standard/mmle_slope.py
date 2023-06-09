# %%
import proplot as pplt
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

# %%
def read_extreme_counts(plev = 50000):
    models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
    # load data
    # different models
    ds = {}
    for model in models:
        ds[model] = xr.open_dataset(
            f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/plev_{plev}_extre_counts_tsurf.nc"
        ).squeeze()

    # random sampled models
    dsr = {}
    for ens_size in np.arange(20, 101, 10):
        dsr[ens_size] = xr.open_dataset(
            f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/plev_{plev}_extreme_counts_tsurf_{str(ens_size)}.nc"
        ).squeeze()

    return ds,dsr

# %%
# calcualte slope
def calc_slope(ds, extr_type, mode):

    x = ds.tsurf.values
    y = ds.extreme_counts.sel(extr_type=extr_type, mode=mode, confidence="true").values

    model = sm.OLS(y, sm.add_constant(x)).fit()
    slope = model.params[1]
    conf_int = model.conf_int()[1]
    return slope, conf_int


# %%
def plot_scatter(ds, models, axs, colors, ensemble_size=None,alpha = 0.7):
    extr_types = ["pos", "neg"]
    modes = ["NAO", "EA"]
    if len(colors) != len(models):
        colors = [colors] * len(models)

    for i, extr_type in enumerate(extr_types):
        ax = axs[i]
        ax.format(title=extr_type)
        for j, model in enumerate(models):
            slope_NAO, conf_int_NAO = calc_slope(ds[model], extr_type, "NAO")
            slope_EA, conf_int_EA = calc_slope(ds[model], extr_type, "EA")

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


#%%
# Create a scatter plot of the slopes for each model and extreme type
def plot_slope(plev = 50000):
    models = ["MPI_GE_onepct","MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
    ds,dsr = read_extreme_counts(plev = plev)
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
    colors = cmap(np.linspace(0.8, 0.3, 9))
    ens_size = np.arange(20, 101, 10)

    plot_scatter(ds, models, axs[0, :], colors_model, ensemble_size=ensemble_size)
    plot_scatter(dsr, np.arange(20, 101, 10), axs[1, :], "tab:red", ensemble_size=ens_size,alpha=0.5)
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

    # plt.savefig(f'/work/mh0033/m300883/Tel_MMLE/docs/source/plots/MMLE/plev_{plev}_slope.png')



# %%
plot_slope(plev=50000)

# %%
plot_slope(plev=30000)
