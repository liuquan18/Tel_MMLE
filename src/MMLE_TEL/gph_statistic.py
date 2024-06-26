# %%
import xarray as xr
import numpy as np
import os
import sys

import src.MMLE_TEL.index_generator as index_generator
import src.Teleconnection.tools as tools
from scipy.stats import linregress


# %%
def read_gph_data(model, remove_ens_mean=True):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/"

    # read gph data
    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + "zg_" + month + "/"
        data_JJA.append(
            index_generator.read_data(
                zg_path, plev=50000, remove_ens_mean=remove_ens_mean
            )
        )
    data = xr.concat(data_JJA, dim="time").sortby("time")

    # drop the year 2100 since it is not complete as a 10-year mean
    data = data.sel(time=slice(None, "2099"))

    data = change_lon_to_180(data)
    return data


# function to change the lontitude from 0-360 to -180-180
def change_lon_to_180(zg):
    zg = zg.assign_coords(lon=(((zg.lon + 180) % 360) - 180)).sortby("lon")
    return zg


# %%


def standardize_arr(arr, standard="first10", dim="com",divide_std=True):
    """
    standard when count the extreme cases.
    """

    if standard == "first10":
        years = np.unique(arr.sortby("time").time.dt.year)
        arr_ref = arr.sel(time=slice(str(years[0]), str(years[9])))
        arr_ref = arr_ref.stack(com=("time", "ens"))
        arr_ref_mean = arr_ref.mean(dim=dim)
        arr_ref_std = arr_ref.std(dim=dim)
        if divide_std:
            arr_standard = (arr - arr_ref_mean) / arr_ref_std
        else:
            arr_standard = arr - arr_ref_mean
    elif standard == "all":
        arr = arr.stack(com=("time", "ens"))
        arr_box_mean = arr.mean(dim=dim)
        arr_box_std = arr.std(dim=dim)
        if divide_std:
            arr_standard = (arr - arr_box_mean) / arr_box_std
        else:
            arr_standard = arr - arr_box_mean
        arr_standard = arr_standard.unstack()
    return arr_standard


def box_arr_cov(box_mean, arr):
    """
    calculate the covariability between a temporal-ens box mean (time*ens x 1)
    and the spatial-temporalens array (time*ens x lat x lon)
    """
    x = box_mean.stack(com=("time", "ens"))
    y = arr.stack(com=("time", "ens"))
    y = y.transpose("com", "lat", "lon").values
    cov = np.cov(x, y.reshape(y.shape[0], -1), rowvar=False)

    # extract covariability between box mean and each grid point
    cov_xy = cov[0, 1:].reshape(y.shape[1], y.shape[2])

    cov_xy = xr.DataArray(
        cov_xy, dims=("lat", "lon"), coords={"lat": arr.lat, "lon": arr.lon}
    )

    return cov_xy

def box_arr_cov_std(box_mean, arr):
    """
    same as box_arr_cov, but set the mean of x and y as zero when calculate the covariance.
    the box_mean and arr should be removed the mean of time and ens of the first10 years.
    """
    x = box_mean.stack(com=("time", "ens"))
    y = arr.stack(com=("time", "ens"))
    y = y.transpose("com", "lat", "lon")

    stadand_cov = xr.apply_ufunc(
        std_cov,
        x,
        y,
        input_core_dims=[["com"], ["com"]],
        output_core_dims=[[]],
        vectorize=True,
    )

    return stadand_cov

def std_cov(x,y):
    """
    calculate covariance between x and y, but do not remove the mean of x and y.
    """
    n = len(x)
    cov = np.sum(x * y) / (n - 1)
    return cov




def stats_arr(arr, statis="std", **kwargs):
    if statis == "std":
        arr = arr.stack(com=("time", "ens"))
        arr_stat = arr.std(dim="com")
    elif statis == "mean":
        arr = arr.stack(com=("time", "ens"))
        arr_stat = arr.mean(dim="com")
    elif statis == "count_pos":
        threshold = kwargs.get("threshold", 1.5)
        arr_com = arr.stack(com=("time", "ens"))
        # count the number of events above 1.5
        arr_stat = arr_com.where(arr_com > threshold).count(dim="com")
    elif statis == "count_neg":
        threshold = kwargs.get("threshold", 1.5)
        arr_com = arr.stack(com=("time", "ens"))
        # count the number of events above 1.5
        arr_stat = arr_com.where(arr_com < -1 * threshold).count(dim="com")
    elif statis == "cov_pos":
        std_cov = kwargs.get("std_cov", False)
        box_mean = box_spatial_mean(arr, 45, 60, -30, 0)
        if std_cov:
            arr_stat = box_arr_cov_std(box_mean, arr)
        else:
            arr_stat = box_arr_cov(box_mean, arr)
    elif statis == "cov_neg":
        std_cov = kwargs.get("std_cov", False)
        box_mean = box_spatial_mean(arr, 60, 75, -75, -50)
        if std_cov:
            arr_stat = box_arr_cov_std(box_mean, arr)
        else:
            arr_stat = box_arr_cov(box_mean, arr)

    else:
        print("please input the correct statistic method")

    return arr_stat


# %%
def box_spatial_mean(xarr, blat, tlat, llon, rlon):
    """
    calculate the spatial mean of the box.
    check if the lat and lon of the eof is from lower to higher.
    input the order of expected box in # bottom lat, top lat, left lon, right lon for the box
    return the order of the box in the eof.
    """
    if xarr.lat[0] > xarr.lat[-1]:  # descending
        start_lat = tlat
        end_lat = blat
    else:  # ascending
        start_lat = blat
        end_lat = tlat

    if xarr.lon[0] > xarr.lon[-1]:
        start_lon = rlon
        end_lon = llon

    else:
        start_lon = llon
        end_lon = rlon

    arr_box = xarr.sel(lat=slice(start_lat, end_lat), lon=slice(start_lon, end_lon))
    weights = np.cos(np.deg2rad(arr_box.lat))
    region_box_weighted = arr_box.weighted(weights)
    arr_box_spatial_mean = region_box_weighted.mean(("lon", "lat"))

    return arr_box_spatial_mean


def box_variability(model):
    zg = read_gph_data(model)

    pos_std = zg.resample(time="10AS-JUN").apply(NAO_pos_center, statis="std")
    neg_std = zg.resample(time="10AS-JUN").apply(NAO_neg_center, statis="std")

    return pos_std, neg_std


def box_mean(model):
    zg = read_gph_data(model)

    pos_mean = zg.resample(time="10AS-JUN").apply(NAO_pos_center, statis="mean")
    neg_mean = zg.resample(time="10AS-JUN").apply(NAO_neg_center, statis="mean")
    return pos_mean, neg_mean


def NAO_pos_center(zg, statis="std"):
    pos_box_spa_mean = box_spatial_mean(zg, 45, 60, -30, 0)
    pos_box_std = stats_arr(pos_box_spa_mean, statis=statis)
    return pos_box_std


def NAO_neg_center(zg, statis="std"):
    neg_box = box_spatial_mean(zg, 60, 75, -75, -50)
    neg_box_std = stats_arr(neg_box, statis=statis)
    return neg_box_std
    # spaital mean with weights


# %%
# define function to calculate slope and p-value
def linregress_ufunc(x):
    slope, _, r_value, p_value, _ = linregress(np.arange(len(x)), x)
    return slope, p_value


# %%
def ens_std(model):
    zg = read_gph_data(model)
    zg.load()
    zg_std = zg.resample(time="10AS-JUN").apply(stats_arr, statis="std")
    return zg_std


def slope_ens_std(model):
    zg_std = ens_std(model)

    result = xr.apply_ufunc(
        linregress_ufunc,
        zg_std,
        input_core_dims=[["time"]],
        output_core_dims=[[], []],
        vectorize=True,
    )

    # convert result to DataArray
    slope_da = result[0].rename("slope")
    pvalue_da = result[1].rename("pvalue")

    # create a new dataset and add the slope_da and pvalue_da as variables
    result_ds = xr.Dataset({"slope": slope_da, "pvalue": pvalue_da})
    return result_ds


# %%


def slope_ens_mean(model):
    zg = read_gph_data(model, remove_ens_mean=False)
    zg.load()

    # calculate the decadal mean of zg to corresponds to the slope of ens std.
    zg_yearly = zg.resample(time="10AS-JUN").mean(dim=("time", "ens"))  # ens mean

    # calculate the slope and pvalue
    result = xr.apply_ufunc(
        linregress_ufunc,
        zg_yearly,
        input_core_dims=[["time"]],
        output_core_dims=[[], []],
        vectorize=True,
    )

    # convert result to DataArray
    slope_da = result[0].rename("slope")
    pvalue_da = result[1].rename("pvalue")

    # create a new dataset and add the slope_da and pvalue_da as variables
    result_ds = xr.Dataset({"slope": slope_da, "pvalue": pvalue_da})
    return result_ds


# %%
def gph_extrc(model, **kwargs):
    zg = read_gph_data(model)
    zg.load()

    # standardize the zg to count the extreme cases (1.5 std)
    standard = kwargs.get("standard", "first10")
    zg_standard = standardize_arr(zg, standard=standard)
    # count the extreme cases every ten years
    zg_count_pos = zg_standard.resample(time="10AS-JUN").apply(
        stats_arr, statis="count_pos"
    )

    zg_count_neg = zg_standard.resample(time="10AS-JUN").apply(
        stats_arr, statis="count_neg"
    )

    return zg_count_pos, zg_count_neg


# %%
def slope_gph_extrc(model, **kwargs):
    zg_count_pos, zg_count_neg = gph_extrc(model, **kwargs)

    # calculate the slope and pvalue
    result_pos = xr.apply_ufunc(
        linregress_ufunc,
        zg_count_pos,
        input_core_dims=[["time"]],
        output_core_dims=[[], []],
        vectorize=True,
    )

    result_neg = xr.apply_ufunc(
        linregress_ufunc,
        zg_count_neg,
        input_core_dims=[["time"]],
        output_core_dims=[[], []],
        vectorize=True,
    )

    # convert result to DataArray
    slope_pos_da = result_pos[0].rename("slope")
    pvalue_pos_da = result_pos[1].rename("pvalue")

    slope_neg_da = result_neg[0].rename("slope")
    pvalue_neg_da = result_neg[1].rename("pvalue")

    # create a new dataset and add the slope_da and pvalue_da as variables
    result_ds = xr.Dataset(
        {
            "slope_pos": slope_pos_da,
            "pvalue_pos": pvalue_pos_da,
            "slope_neg": slope_neg_da,
            "pvalue_neg": pvalue_neg_da,
        }
    )
    return result_ds


# %%
def box_cov(model, **kwargs):
    zg = read_gph_data(model)
    zg.load()

    pos_cov = zg.resample(time="10AS-JUN").apply(stats_arr, statis="cov_pos")
    neg_cov = zg.resample(time="10AS-JUN").apply(stats_arr, statis="cov_neg")

    return pos_cov, neg_cov


# %%
def slope_box_cov(model, **kwargs):
    pos_cov, neg_cov = box_cov(model, **kwargs)
    # calculate the slope and pvalue
    result_pos = xr.apply_ufunc(
        linregress_ufunc,
        pos_cov,
        input_core_dims=[["time"]],
        output_core_dims=[[], []],
        vectorize=True,
    )

    result_neg = xr.apply_ufunc(
        linregress_ufunc,
        neg_cov,
        input_core_dims=[["time"]],
        output_core_dims=[[], []],
        vectorize=True,
    )

    # convert result to DataArray
    slope_pos_da = result_pos[0].rename("slope")
    pvalue_pos_da = result_pos[1].rename("pvalue")

    slope_neg_da = result_neg[0].rename("slope")
    pvalue_neg_da = result_neg[1].rename("pvalue")

    # create a new dataset and add the slope_da and pvalue_da as variables
    result_ds = xr.Dataset(
        {
            "slope_pos": slope_pos_da,
            "pvalue_pos": pvalue_pos_da,
            "slope_neg": slope_neg_da,
            "pvalue_neg": pvalue_neg_da,
        }
    )
    return result_ds
