# %%
import importlib
import matplotlib.pyplot as plt
import xarray as xr
import os
import src.extreme.extreme_ci as extreme
import src.composite.composite as composite
import numpy as np
import glob

# %%
import importlib

# %%


# %%
def extreme_counts(
    model,
    fixed_pattern="decade",
    standard="first",
    tsurf="ens_fld_year_mean",
    plev=50000,
    season="JJAS",
):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/"
    prefix = f"plev_{plev}_{fixed_pattern}_{standard}_{season}_"
    eof_dir = odir + "EOF_result/" + prefix + "eof_result.nc"

    # to dir
    extr_counts_dir = odir + "extreme_count/" + prefix + "extre_counts.nc"

    # eof
    eof_result = xr.open_dataset(eof_dir)
    # extreme counts
    # check if the extr_counts_dir exists
    extrc = extreme.decadal_extrc(eof_result.pc, ci="bootstrap")

    # save the result
    try:
        extrc.to_netcdf(extr_counts_dir)
    except PermissionError:
        # remove the file
        os.remove(extr_counts_dir)
        extrc.to_netcdf(extr_counts_dir)


# %%
def decadal_tsurf(
    model,
    tsurf="ens_fld_year_mean",
):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/"
    tsurf_dir = odir + "ts_processed/" + tsurf + ".nc"

    # to dir
    t_surf_mean_dir = odir + "extreme_count/" + tsurf + ".nc"

    # tsurf
    temperature = read_spmean_tsurf(tsurf_dir)

    # tsurf mean
    tsurf_mean = extreme.decade_tsurf(temperature)

    # save the result
    try:
        tsurf_mean.to_netcdf(t_surf_mean_dir)
    except PermissionError:
        os.remove(t_surf_mean_dir)
        tsurf_mean.to_netcdf(t_surf_mean_dir)


# %%
# read tsurf
def read_spmean_tsurf(tsurf_dir):
    tsurf = xr.open_dataset(tsurf_dir)

    # read the array either tusrf, ts, or tas
    try:
        tsurf_arr = tsurf.tsurf.squeeze()
    except AttributeError:
        tsurf_arr = tsurf.ts.squeeze()
    # change the temp time into datetime64
    try:
        tsurf_arr["time"] = tsurf_arr.indexes["time"].to_datetimeindex()
    except AttributeError:
        pass

    return tsurf_arr


def read_var_data(field_var_dir, remove_ensmean=True):
    print(" reading the field data...")

    all_ens_lists = sorted(
        glob.glob(field_var_dir + "*.nc")
    )  # to make sure that the order of ensemble members is fixed
    var_data = xr.open_mfdataset(all_ens_lists, combine="nested", concat_dim="ens")
    var_data["ens"] = np.arange(var_data.ens.size)

    try:
        var_data["time"] = var_data["time"].astype("datetime64[ns]")
    except TypeError:
        pass

    # demean the ensemble mean
    if remove_ensmean:
        var_data = var_data - var_data.mean(dim="ens")
    return var_data


def split_first_last(eof_result):
    times = eof_result.time
    years = np.unique(times.dt.year)
    first_years = years[:10]
    last_years = years[-10:]

    eof_first = eof_result.isel(decade=0).sel(
        time=eof_result["time.year"].isin(first_years)
    )
    eof_last = eof_result.isel(decade=-1).sel(
        time=eof_result["time.year"].isin(last_years)
    )
    return eof_first, eof_last


# %%
def extreme_counts_profile(model, standard="first", season="MJJA"):
    eof_dir = (
        f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"
        + "troposphere_ind_decade_"
        + standard
        + "_"
        + season
        + "_eof_result.nc"
    )

    eof_result = xr.open_dataset(eof_dir)

    # PCS of the first and last decade
    first_pc = eof_result.pc.isel(time=slice(0, 10))
    last_pc = eof_result.pc.isel(time=slice(-10, None))

    if first_pc.time.size != 10 or last_pc.time.size != 10:
        raise ValueError("the time size is not 10")

    print("calculating the extreme event count")
    # extreme counts of first and last 10 decades
    first_count = extreme.extreme_count_xr(first_pc, ci="bootstrap")
    last_count = extreme.extreme_count_xr(last_pc, ci="bootstrap")

    # save the result
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    prefix = f"troposphere_ind_decade_{standard}_{season}_"
    first_count.to_netcdf(odir + prefix + "first_count.nc")
    last_count.to_netcdf(odir + prefix + "last_count.nc")


# %%
# composite analysis of surface temperature in terms of different extreme events
def composite_analysis(
    model,
    fixed_pattern="decade",
    standard="first",
    plev=50000,
    index_season="MJJA",
    var_season="MJJA",
    var_data=None,
    reduction="mean",  # mean_same_number, mean_weighted, mean
    threshold=1.5,
    var_name="ts",
    **kwargs,
):
    """
    tfield can be 'same' or 'next'
    """
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/"
    prefix = f"plev_{plev}_{fixed_pattern}_{standard}_{index_season}_"
    eof_dir = odir + "EOF_result/" + prefix + "eof_result.nc"
    field_var_dir = odir + f"{var_name}_{var_season}/"

    # to dir
    composite_dir = f"{odir}composite/{prefix}{var_season}_first_last_{var_name}_composite_{reduction}.nc"

    if var_data is None:
        var_data = read_var_data(field_var_dir)
        if var_name == "ts":
            try:
                var_data = var_data["tsurf"]
            except KeyError:
                var_data = var_data["ts"]
        elif var_name == "pr":
            try:
                var_data = var_data["pr"]
            except KeyError:
                var_data = var_data["precip"]
    else:
        pass

    var_data.load()

    # eof
    eof_result = xr.open_dataset(eof_dir)
    eof_first, eof_last = split_first_last(eof_result)

    # select the first and last 10 decades

    first_index = eof_first.pc
    last_index = eof_last.pc

    print(f" compositing {reduction} of the {var_name} data...")
    composite_mean = composite.first_last_extreme_composite(
        first_index,
        last_index,
        var_data,
        threshold=threshold,
        reduction=reduction,
        return_diff=True,
        **kwargs,
    )
    # save the result
    try:
        composite_mean.to_netcdf(composite_dir)
    except PermissionError:
        os.remove(composite_dir)
        composite_mean.to_netcdf(composite_dir)
