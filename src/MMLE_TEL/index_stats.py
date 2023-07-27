# %%
import src.MMLE_TEL.story_line as story_line
import importlib
import matplotlib.pyplot as plt
import xarray as xr
import os
import src.extreme.extreme_ci as extreme
import src.composite.field_composite as composite


#%%
import importlib

importlib.reload(story_line)

# %%

#%%
def extreme_counts_tsurf(
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
    tsurf_dir = odir + "ts_processed/" + tsurf + ".nc"

    # to dir
    extr_counts_dir = odir + "extreme_count/" + prefix + "extre_counts.nc"
    t_surf_mean_dir = odir + "extreme_count/" + tsurf + ".nc"

    # eof
    eof_result = xr.open_dataset(eof_dir)
    # tsurf
    temperature = read_tsurf(tsurf_dir)

    # extreme counts
    # check if the extr_counts_dir exists
    try:
        extrc = xr.open_dataset(extr_counts_dir).pc

    except FileNotFoundError:
        extrc = extreme.decadal_extrc(eof_result.pc, ci="bootstrap")

    # tsurf mean
    tsurf_mean = extreme.decade_tsurf(extrc, temperature)

    # save the result
    try:
        extrc.to_netcdf(extr_counts_dir)
    except PermissionError:
        pass
    tsurf_mean.to_netcdf(t_surf_mean_dir)

    return extrc, tsurf_mean


#%%
# read tsurf
def read_tsurf(tsurf_dir):
    tsurf = xr.open_dataset(tsurf_dir)

    # read the array either tusrf, ts, or tas
    try:
        tsurf_arr = tsurf.tsurf.squeeze()
    except AttributeError:
        try:
            tsurf_arr = tsurf.ts.squeeze()
        except AttributeError:
            tsurf_arr = tsurf.tas.squeeze()
    # change the temp time into datetime64
    try:
        tsurf_arr["time"] = tsurf_arr.indexes["time"].to_datetimeindex()
    except AttributeError:
        pass

    return tsurf_arr


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

#%%
# composite analysis of surface temperature in terms of different extreme events
def composite_analysis(
        model,
        fixed_pattern="decade",
        standard="first",
        plev=50000,
        index_season="MJJA",
        tsurf_season="MJJA",
        reduction="mean",
        threshold=1.5,):
    """
    tfield can be 'same' or 'next'
    """
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/"
    prefix = f"plev_{plev}_{fixed_pattern}_{standard}_{index_season}_"
    eof_dir = odir + "EOF_result/" + prefix + "eof_result.nc"
    field_tsurf_dir = odir + f"ts_{tsurf_season}/"

    # to dir
    first_composite_dir = f"{odir}composite/{prefix}{tsurf_season}_first_composite.nc"
    last_composite_dir = f"{odir}composite/{prefix}{tsurf_season}_last_composite.nc"

    print(" reading the tsurf field data...")
    try:
        var_data = xr.open_dataset(field_tsurf_dir + "all_ens_tsurf.nc")
    except FileNotFoundError:
        var_data = xr.open_mfdataset(field_tsurf_dir + "*.nc",combine='nested',concat_dim='ens')

    try:
        var_data = var_data.tsurf
    except AttributeError:
        try:
            var_data = var_data.ts
        except AttributeError:
            var_data = var_data.tas
    try:
        var_data['time'] = var_data['time'].astype('datetime64[ns]')
    except TypeError:
        pass
    var_data = var_data - var_data.mean(dim="ens")


    # eof
    eof_result = xr.open_dataset(eof_dir)
    index = eof_result.pc

    # select the first and last 10 decades
    first_index = index.isel(time=slice(0, 10))
    last_index = index.isel(time=slice(-10, None))

    print(" compositing the tsurf data...")
    first_var = composite.Tel_field_composite(
        first_index, var_data, threshold=threshold, reduction=reduction
    )
    last_var = composite.Tel_field_composite(
        last_index, var_data, threshold=threshold, reduction=reduction
    )


    # save the result
    first_var.to_netcdf(first_composite_dir)
    last_var.to_netcdf(last_composite_dir)