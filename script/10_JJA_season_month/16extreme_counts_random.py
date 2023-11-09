# %%
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf
import src.extreme.extreme_ci as extreme
import xarray as xr
import os

# %%
import importlib


# %%
def read_eof(
    ens_size,
    model="MPI_GE",
    plev=50000,
    standard="first",
    tsurf="ens_fld_year_mean",
    season="JJA",
    fixedPattern="decade",
):
    print(f"read data of {model} {tsurf} {season} {standard} {str(ens_size)}")

    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/EOF_result/"
    standard_name = (
        "plev_"
        + str(plev)
        + "_"
        + fixedPattern
        + "_"
        + season
        + "_"
        + standard
        + "_"
        + str(ens_size)
        + "_eof_result.nc"
    )
    eof_dir = odir + standard_name
    eof_result = xr.open_dataset(eof_dir)

    return eof_result


# %%
# count the extreme event and tsurf
def extreme_counts(eof_result):
    extrc = extreme.decadal_extrc(eof_result.pc, ci=None)
    return extrc


# %%
def process_data(ens_size, model="MPI_GE"):
    plev = 50000
    fixed_pattern = "decade"
    eof_result = read_eof(ens_size, model=model)
    season = "JJA"
    standard = "first"

    extr_counts = extreme_counts(eof_result)
    try:
        extr_counts.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/extreme_count/plev_{plev}_{fixed_pattern}_{season}_{standard}_{str(ens_size)}_extre_counts.nc"
        )
    except PermissionError:
        os.remove(
            f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/extreme_count/plev_{plev}_{fixed_pattern}_{season}_{standard}_{str(ens_size)}_extre_counts.nc"
        )
        extr_counts.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/extreme_count/plev_{plev}_{fixed_pattern}_{season}_{standard}_{str(ens_size)}_extre_counts.nc"
        )


# %%

ens_sizes = [20, 30, 40, 50]
for ens_size in ens_sizes:
    process_data(ens_size, model="MPI_GE")

# %%
