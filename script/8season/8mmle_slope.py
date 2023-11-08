# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

# %%
import src.plots.extreme_plot as ext_plot
import importlib

importlib.reload(ext_plot)


# %%
def read_extreme_counts(**kwargs):
    plev = kwargs.get("plev", 50000)
    standard = kwargs.get("standard", "first")
    tsurf = kwargs.get("tsurf", "ens_fld_year_mean")
    fixed_pattern = kwargs.get("fixed_pattern", "decade")
    season = kwargs.get("season", "JJAS")

    models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
    # load data
    # different models
    extrs = {}
    tsurfs = {}
    for model in models:
        odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
        prefix = f"plev_{plev}_{fixed_pattern}_{standard}_{season}_"
        extrs[model] = xr.open_dataset(f"{odir}{prefix}extre_counts.nc").squeeze()
        tsurfs[model] = xr.open_dataset(f"{odir}{tsurf}.nc")
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
            odir
            + f"plev_{plev}_{fixed_pattern}_{season}_{standard}_{str(ens_size)}_extre_counts.nc"
        ).squeeze()
        tsurfs_rand[ens_size] = xr.open_dataset(
            f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/{tsurf}.nc"
        ).tsurf

    return extrs, tsurfs, extrs_rand, tsurfs_rand


# %%
# Create a scatter plot of the slopes for each model and extreme type
def mmle_slope(**kwargs):
    tsurf = kwargs.get("tsurf", "ens_fld_year_mean")
    time = kwargs.get("time", "all")
    season = kwargs.get("season", "JJAS")
    plev = kwargs.get("plev", 50000)

    extrs, tsurfs, extrs_rand, tsurfs_rand = read_extreme_counts(**kwargs)
    ext_plot.mmle_slope_scatter(
        extrs, tsurfs, extrs_rand, tsurfs_rand, tsurf=tsurf, time=time
    )
    # plt.savefig(
    #     f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Winter_Summer/plev_{plev}_extrc_{tsurf}_{time}_{season}_slope.png"
    # )


# %%
def slope_diff_tsurf(**kwargs):
    plev = kwargs.get("plev", 50000)
    season = kwargs.get("season", "JJAS")
    time = kwargs.get("time", "all")

    extrs, tsurf_gmst, _, _ = read_extreme_counts(**kwargs)
    extrs, NA_tsurf, _, _ = read_extreme_counts(**kwargs)
    extrs, tropical_arctic_gradient, _, _ = read_extreme_counts(**kwargs)

    ext_plot.slope_diff_tsurf(
        extrs, tsurf_gmst, NA_tsurf, tropical_arctic_gradient, time=time
    )

    plt.savefig(
        f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Winter_Summer/plev_{plev}_extrc_tsurf_all_{time}_slope_{season}_.png"
    )


# %%
mmle_slope(
    plev=50000,
    standard="first",
    tsurf="ens_fld_year_mean",
    time="1960-01-01",
    season="DJFM",
)

# %%
mmle_slope(
    plev=50000,
    standard="first",
    tsurf="ens_fld_year_mean",
    time="1960-01-01",
    season="JJAS",
)

# %%
mmle_slope(
    plev=50000,
    standard="first",
    tsurf="ens_fld_year_mean",
    time="1960-01-01",
    season="MAM",
)
# %%


# %%
extrs, tsurfs, extrs_rand, tsurfs_rand = read_extreme_counts(
    plev=50000, standard="first", tsurf="ens_fld_year_mean", season="JJAS"
)
# %%
importlib.reload(ext_plot)
fig = ext_plot.mmle_line_plot(
    extrs, tsurfs, extrs_rand, tsurfs_rand, tsurf="ens_fld_year_mean", time="1960-01-01"
)
# fig.savefig('/work/mh0033/m300883/Tel_MMLE/docs/source/plots/slides_IUGG/mmle_line_plot_GFDL_CM3.png', facecolor=fig.get_facecolor(), edgecolor='none')
# %%
