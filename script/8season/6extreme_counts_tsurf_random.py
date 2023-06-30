# %%
import src.MMLE_TEL.story_line as story_line
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf
import src.extreme.extreme_ci as extreme
import xarray as xr

#%%
import importlib

importlib.reload(story_line)


# %%
def read_data(
    ens_size,
    model="MPI_GE_onepct",
    plev=50000,
    standard="first",
    tsurf="ens_fld_year_mean",
    season = 'MAM',
    fixedPattern = "decade",
):
    eof_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/EOF_result/plev_{plev}_{fixedPattern}_{season}_{standard}_{str(ens_size)}_eof_result.nc"
    eof_result = xr.open_dataset(eof_dir)

    tsurf_dir = (
        f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/ts_processed/{tsurf}.nc"
    )
    tsurf = xr.open_dataset(tsurf_dir).tsurf

    extrc_dir = f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/plev_{plev}_{fixedPattern}_{season}_{standard}_{str(ens_size)}_extre_counts.nc"
    try:
        extrc = xr.open_dataset(extrc_dir).pc
    except FileNotFoundError:
        extrc = None
    return eof_result, tsurf, extrc


# %%
# count the extreme event and tsurf
def extreme_counts(eof_result, tsurf, extrc=None):

    tsurf_increase = tsurf - tsurf[0]
    if extrc is None:
        ext_counts, t_surf_mean = extreme.decadal_extrc_tsurf(
            eof_result.pc, temp=tsurf_increase, ci="bootstrap"
        )
    else:
        ext_counts, t_surf_mean = extreme.decadal_extrc_tsurf(
            eof_result.pc, ext_counts_xr=extrc, temp=tsurf_increase, ci="bootstrap"
        )
    return ext_counts, t_surf_mean


# %%
from multiprocessing import Pool

ens_sizes = np.arange(20, 101, 10)
tsurfs = ["ens_fld_year_mean"]#, "NA_tsurf", "tropical_arctic_gradient"]
standards = ["first"]#, "temporal_ens"]


def process_data(ens_size, tsurf, standard,season):
    plev = 50000
    fixed_pattern = "decade"
    season = "MAM"
    eof_result, temperature, extrc = read_data(ens_size, tsurf=tsurf, season = 'MAM')
    # close the files red
    temperature.close()

    extr_counts, tsurf_mean = extreme_counts(eof_result, temperature, extrc=extrc)
    try:
        extr_counts.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/plev_{plev}_{fixed_pattern}_{season}_{standard}_{str(ens_size)}_extre_counts.nc"
        )
    except PermissionError:
        pass
    try:
        tsurf_mean.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/{tsurf}.nc"
        )
    except PermissionError:
        pass


if __name__ == "__main__":
    with Pool() as p:
        p.starmap(
            process_data,
            [
                (ens_size, tsurf, standard,season)
                for ens_size in ens_sizes
                for tsurf in tsurfs
                for standard in standards
                for season in ['MAM']
            ],
        )


# %%
