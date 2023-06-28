# %%
import proplot as pplt
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
#%%
import src.plots.extreme_plot as ext_plot


# %%
def read_extreme_counts(plev = 50000,standard = 'first',tsurf = 'ens_fld_year_mean',fixed_pattern = 'decade',season = 'JJAS'):
    models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
    # load data
    # different models
    extrs = {}
    tsurfs = {}
    for model in models:
        odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
        prefix = f"plev_{plev}_{fixed_pattern}_{standard}_{season}_"
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


#%%
# Create a scatter plot of the slopes for each model and extreme type
def plot_slope(plev = 50000,standard = 'first',tsurf = 'ens_fld_year_mean',time = 'all'):
    ext_plot.mmle_slope_scatter(plev = plev,standard = standard,tsurf = tsurf,time = time)
    plt.savefig(f'/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Winter_Summer/plev_{plev}_extrc_{tsurf}_{time}_slope.png')

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

    axs[2,:].format(
        # reverse the x and y axis
        xlim = (1.9,-3),
        ylim = (1.9,-3.5)
    )


    plt.savefig(f'/work/mh0033/m300883/Tel_MMLE/docs/source/plots/MMLE/plev_{plev}_extrc_tsurf_all_{time}_slope.png')

#%%




# %%
plot_slope(plev=50000,tsurf='ens_fld_year_ocean_mean',standard='first')
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
# write the slope plots into a md file
def create_doc():
    print("creating the doc")
    with open(
        "/work/mh0033/m300883/Tel_MMLE/docs/source/MMLE_slope.md", "w"
    ) as f:
        f.write("# Slopes of extreme counts vs. temperature\n")

        f.write("## GMST with all the time\n")
        f.write("![GMST with all the time](plots/MMLE/plev_50000_extrc_ens_fld_year_mean_all_slope.png)\n")

        f.write("## GMST with time after 1950\n")
        f.write("![GMST with time after 1950](plots/MMLE/plev_50000_extrc_ens_fld_year_mean_1950-01-01_slope.png)\n")

        f.write("## all the time with different tsurf\n")
        f.write("![all the time with different tsurf](plots/MMLE/plev_50000_extrc_tsurf_all_all_slope.png)\n")

        f.write("## time after 1950 with different tsurf\n")
        f.write("![time after 1950 with different tsurf](plots/MMLE/plev_50000_extrc_tsurf_all_1950-01-01_slope.png)\n")
# %%
