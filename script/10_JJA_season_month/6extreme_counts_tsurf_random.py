# %%
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf
import src.extreme.extreme_ci as extreme
import xarray as xr
import os

#%%
import importlib



# %%
def read_data(
    ens_size,
    model="MPI_GE_onepct",
    plev=50000,
    standard="first",
    tsurf="ens_fld_year_mean",
    season = 'JJA',
    fixedPattern = "decade_mpi",
):
    print(f"read data of {model} {tsurf} {season} {standard} {str(ens_size)}")
    eof_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/EOF_result/plev_{plev}_{fixedPattern}_{season}_{standard}_{str(ens_size)}_eof_result.nc"
    eof_result = xr.open_dataset(eof_dir)
    return eof_result


# %%
# count the extreme event and tsurf
def extreme_counts(eof_result):
    ext_counts = extreme.decadal_extrc(
        eof_result.pc,  ci="bootstrap"
    )
    return ext_counts

def process_data(ens_size, tsurf, standard,season):
    plev = 50000
    fixed_pattern = "decade_mpi"
    season = season
    eof_result= read_data(ens_size, tsurf=tsurf, season = season)
    # close the files red

    extr_counts = extreme_counts(eof_result)
    try:
        extr_counts.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/plev_{plev}_{fixed_pattern}_{season}_{standard}_{str(ens_size)}_extre_counts.nc"
        )
    except PermissionError:
        os.remove(
            f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/plev_{plev}_{fixed_pattern}_{season}_{standard}_{str(ens_size)}_extre_counts.nc"
        )
        extr_counts.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/plev_{plev}_{fixed_pattern}_{season}_{standard}_{str(ens_size)}_extre_counts.nc"
        )
#%%

if __name__ == "__main__":
    from multiprocessing import Pool

    ens_sizes = np.arange(20, 101, 10)
    tsurfs = ["ens_fld_year_mean"]#, "NA_tsurf", "tropical_arctic_gradient"]
    standards = ["first"]#, "temporal_ens"]


    with Pool() as p:
        p.starmap(
            process_data,
            [
                (ens_size, tsurf, standard,season)
                for ens_size in ens_sizes
                for tsurf in tsurfs
                for standard in standards
                for season in ['JJA']
            ],
        )


# %%
