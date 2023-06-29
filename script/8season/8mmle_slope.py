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
def mmle_slope(plev = 50000,standard = 'first',tsurf = 'ens_fld_year_mean',time = 'all',season = 'JJAS'):
    extrs, tsurfs, extrs_rand, tsurfs_rand  = read_extreme_counts(plev = plev,standard = standard,tsurf = tsurf,season = season)
    ext_plot.mmle_slope_scatter(extrs, tsurfs, extrs_rand, tsurfs_rand,tsurf = tsurf,time = time)
    plt.savefig(f'/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Winter_Summer/plev_{plev}_extrc_{tsurf}_{time}_{season}_slope.png')

#%%
def slope_diff_tsurf(plev = 50000,standard = 'first',tsurf = 'ens_fld_year_mean',time = 'all',season = 'JJAS'):
    extrs,tsurf_gmst,_,_ = read_extreme_counts(plev = plev,standard = standard,tsurf = 'ens_fld_year_mean',season = season)
    extrs,NA_tsurf,_,_ = read_extreme_counts(plev = plev,standard = standard,tsurf = 'NA_tsurf',season = season)
    extrs,tropical_arctic_gradient,_,_ = read_extreme_counts(plev = plev,standard = standard,tsurf = 'tropical_arctic_gradient', season = season)

    ext_plot.slope_diff_tsurf(extrs,tsurf_gmst,NA_tsurf,tropical_arctic_gradient,time = time)

    plt.savefig(f'/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Winter_Summer/plev_{plev}_extrc_tsurf_all_{time}_slope_{season}_.png')

#%%
mmle_slope(plev = 50000,standard = 'first',tsurf = 'ens_fld_year_mean',time = 'all',season = 'MAM')


